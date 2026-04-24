#ifndef COMPUTE_SUBGRAPH_H
#define COMPUTE_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/termination_data.h"
#include "../state/barrier_state.h"
#include "../task/barrier_tasks.h"
#include "../task/change_timestep_task.h"
#include "../task/div_exchange_task.h"
#include "../task/forkable_parallel_compute_lane.h"
#include "../task/parallel_compute_lane.h"
#include "../tool/thread_budget.h"

/// Build the unified Compute sub-graph (non-CC_IBM only).
///
/// Merges predictor and corrector into two parallel lanes:
///   ForkableLane: P1, P2 (fork), C1, C2 (fork), C3 — HEAVY phases + AsyncWorker
///   ParallelLane: P3, P4, P5, C4, C5, C6, C7, C8 — MEDIUM/LIGHT phases
///
/// Fork2 pipeline: when !HT3D, C3 emits PreDivP1 directly to ParallelLane,
/// enabling per-mesh overlap of C3 radiation with C4 DivP1.
/// When HT3D, C3 emits PostWallBC to Exchange(6) barrier → PreDivP1 → C4.
///
/// Template parameter PressureTag:
///   PredictorPressure when pressure subgraph active (P3/C6 exit subgraph),
///   PredDivP2Out otherwise (P3 → internal barrier).
template<MeshState PressureTag = MeshState::PredDivP2Out>
inline auto buildComputeSubgraphImpl(int nmeshes, const ThreadBudget &budget) {
    constexpr bool useParallelPressure = (PressureTag != MeshState::PredDivP2Out);

    auto subgraph = std::make_shared<hh::Graph<9,
        // 9 inputs:
        MeshData<>,                                  // P1 from TimestepTask
        MeshData<MeshState::PostPredExch>,           // from exchange → ExchInsPart barrier → P2
        MeshData<MeshState::PostPredVelExch>,        // from exchange → P5
        MeshData<MeshState::PostCorrStep1>,          // from exchange → C2
        MeshData<MeshState::PostParticleOps>,        // from exchange → HVAC barrier
        MeshData<MeshState::PostRadExch>,            // from exchange → Fork2Join
        MeshData<MeshState::PredictorPressure>,      // from pressure subgraph → P4
        MeshData<MeshState::CorrectorPressure>,      // from pressure subgraph → C7
        TerminationData,
        // 8 outputs:
        MeshData<MeshState::MeshExch1>,              // P1 → exchange
        MeshData<MeshState::MeshExch3>,              // ChangeTimeStep → exchange
        MeshData<MeshState::MeshExch4>,              // C1 → exchange
        MeshData<MeshState::MeshExch7>,              // C2 → exchange
        MeshData<MeshState::MeshExch2>,              // C3 → exchange
        MeshData<>,                                  // C8 → TimestepTask
        MeshData<MeshState::PredictorPressure>,      // P3 → pressure subgraph
        MeshData<MeshState::CorrectorPressure>       // C6 → pressure subgraph
    >>("Compute");

    bool ht3d = fds_is_ht3d() != 0;
    bool useAsyncWorker = budget.useAsyncWorker;

    // --- Two parallel lanes ---
    auto forkableLane = std::make_shared<ForkableParallelComputeLane>(
        budget.forkableLane, ht3d, useAsyncWorker);
    auto parallelLane = std::make_shared<ParallelComputeLane<PressureTag>>(budget.parallelLane);

    // DivExchangeTask: shared between predictor (P2→P3) and corrector (C5→C6).
    // Use Default tag when no pressure subgraph so useParallelPressure_ is correct.
    constexpr MeshState DivPressureTag =
        useParallelPressure ? MeshState::PredictorPressure : MeshState::Default;
    auto divExchangeTask = std::make_shared<DivExchangeTask<DivPressureTag>>(
        nmeshes, budget.divExchange);

    // --- Predictor barriers ---
    auto postPredExchBarrier = makeBarrierSM<MeshState::PostPredExch>(
        nmeshes, "ExchInsPart+Hvac",
        "EXCH_INS_PART\\nHVAC_CALC",
        [](auto& meshes) {
            fds_exchange_inserted_particles();
            fds_hvac_calc(meshes[0]->t, meshes[0]->dt, 1);
        });

    auto changeTimeStepTask = std::make_shared<ChangeTimeStepTask>(
        nmeshes, budget.retryMomDiv, false);

    auto phaseTransTask = std::make_shared<
        PhaseTransitionTask<MeshState::PredFinalOutput, MeshState::CorrInput>>(nmeshes);

    // --- Corrector barriers ---
    auto hvacBarrier = makeTerminableEagerDualMeshBarrier<
        MeshState::PostParticleOps, MeshState::PostHvac>(
        nmeshes, "HvacCalc",
        "HVAC_CALC",
        [](auto& meshes) {
            fds_hvac_calc(meshes[0]->t, meshes[0]->dt, 1);
        });

    // Fork2Join: PostDivP1 (primary from C4) + PostRadExch (secondary from exchange)
    auto fork2JoinTask = std::make_shared<TerminableDualMeshJoinTask<
        MeshState::PostRadExch, MeshState::PostDivJoin, MeshState::PostDivP1>>(
        nmeshes, "Fork2Join");

    int wallIncrement = fds_get_wall_increment();
    auto corrFinalBarrier = makeTerminableRetaggingBarrier<
        MeshState::PostVelCorr, MeshState::PostCorrFinalBarrier>(
        nmeshes, "CorrFinalBarrier",
        "RESET_WALL\\nMESH_EXCHANGE(6)",
        [wallIncrement](auto& meshes) {
            if (meshes[0]->wall_counter == wallIncrement) {
                fds_set_wall_counter(0);
            }
            fds_mesh_exchange(6);
        });

    // =====================================================================
    // PREDICTOR WIRING
    // =====================================================================

    // P1: subgraph input → ForkableLane
    subgraph->template input<MeshData<>>(forkableLane);

    // P1 out: ForkableLane → MeshExch1 → exchange
    subgraph->template output<MeshData<MeshState::MeshExch1>>(forkableLane);

    // PostPredExch → ExchInsPart+Hvac barrier → P2
    subgraph->template input<MeshData<MeshState::PostPredExch>>(postPredExchBarrier);
    subgraph->template edge<MeshData<MeshState::PostPredExch>>(
        postPredExchBarrier, forkableLane);

    // P2 → DivExch → DivExchangeTask → DivPart2 → P3
    subgraph->template edge<MeshData<MeshState::DivExch>>(forkableLane, divExchangeTask);
    subgraph->template edge<MeshData<MeshState::DivPart2>>(divExchangeTask, parallelLane);

    // P3 → Pressure → P4
    if constexpr (useParallelPressure) {
        // P3 exits subgraph → pressure subgraph → re-enters as PredictorPressure → P4
        subgraph->template output<MeshData<MeshState::PredictorPressure>>(parallelLane);
        subgraph->template input<MeshData<MeshState::PredictorPressure>>(parallelLane);
    } else {
        // Internal barrier: PredDivP2Out → pressure iteration → PredDivP2Out
        auto predPressureBarrier = makeTerminableRetaggingBarrier<
            MeshState::PredDivP2Out, MeshState::PredDivP2Out>(
            nmeshes, "PredPressure",
            "PRESSURE_ITERATION\\nINIT_CHANGE_TIME_STEP",
            [](auto& meshes) {
                fds_pressure_iteration(meshes[0]->t, meshes[0]->dt);
                fds_init_change_time_step(meshes[0]->dt);
            });
        subgraph->template edge<MeshData<MeshState::PredDivP2Out>>(
            parallelLane, predPressureBarrier);
        subgraph->template edge<MeshData<MeshState::PredDivP2Out>>(
            predPressureBarrier, parallelLane);
        subgraph->template input<TerminationData>(predPressureBarrier);
    }

    // P4 → PostVelPred → ChangeTimeStep → MeshExch3 → exchange
    subgraph->template edge<MeshData<MeshState::PostVelPred>>(
        parallelLane, changeTimeStepTask);
    subgraph->template output<MeshData<MeshState::MeshExch3>>(changeTimeStepTask);

    // PostPredVelExch → P5
    subgraph->template input<MeshData<MeshState::PostPredVelExch>>(parallelLane);

    // P5 → PredFinalOutput → PhaseTransition → CorrInput → C1
    subgraph->template edge<MeshData<MeshState::PredFinalOutput>>(
        parallelLane, phaseTransTask);
    subgraph->template edge<MeshData<MeshState::CorrInput>>(
        phaseTransTask, forkableLane);

    // =====================================================================
    // CORRECTOR WIRING
    // =====================================================================

    // C1 → MeshExch4 → exchange
    subgraph->template output<MeshData<MeshState::MeshExch4>>(forkableLane);

    // PostCorrStep1 → C2
    subgraph->template input<MeshData<MeshState::PostCorrStep1>>(forkableLane);

    // C2 → Default → HVAC barrier (pre-ParticleOps)
    subgraph->template edge<MeshData<>>(forkableLane, hvacBarrier);
    // C2 → MeshExch7 → exchange (post-ParticleOps)
    subgraph->template output<MeshData<MeshState::MeshExch7>>(forkableLane);
    // PostParticleOps → HVAC barrier (secondary input)
    subgraph->template input<MeshData<MeshState::PostParticleOps>>(hvacBarrier);
    // HVAC barrier → PostHvac → C3
    subgraph->template edge<MeshData<MeshState::PostHvac>>(hvacBarrier, forkableLane);

    // C3 → MeshExch2 → radiation exchange
    subgraph->template output<MeshData<MeshState::MeshExch2>>(forkableLane);

    // C3 → Fork2 path
    if (ht3d) {
        // HT3D: PostWallBC → Exchange(6) barrier → PreDivP1 → C4
        auto exch6Barrier = makeRetaggingBarrier<
            MeshState::PostWallBC, MeshState::PreDivP1>(
            nmeshes, "MeshExch6a",
            "MESH_EXCHANGE(6) [HT3D]",
            [](auto& meshes) {
                if (meshes[0]->call_ht_1d) { fds_mesh_exchange(6); }
            });
        subgraph->template edge<MeshData<MeshState::PostWallBC>>(
            forkableLane, exch6Barrier);
        subgraph->template edge<MeshData<MeshState::PreDivP1>>(
            exch6Barrier, parallelLane);
    } else {
        // !HT3D: ForkableLane emits PreDivP1 directly → per-mesh pipeline
        subgraph->template edge<MeshData<MeshState::PreDivP1>>(
            forkableLane, parallelLane);
    }

    // C4 → PostDivP1 → Fork2Join (primary)
    subgraph->template edge<MeshData<MeshState::PostDivP1>>(
        parallelLane, fork2JoinTask);
    // PostRadExch → Fork2Join (secondary, from exchange)
    subgraph->template input<MeshData<MeshState::PostRadExch>>(fork2JoinTask);

    // Fork2Join → PostDivJoin → C5
    subgraph->template edge<MeshData<MeshState::PostDivJoin>>(
        fork2JoinTask, parallelLane);

    // C5 → DivExch → DivExchangeTask (shared with pred) → DivPart2 → C6
    subgraph->template edge<MeshData<MeshState::DivExch>>(
        parallelLane, divExchangeTask);
    // DivPart2 edge already wired above (divExchangeTask → parallelLane)

    // C6 → Pressure → C7
    if constexpr (useParallelPressure) {
        subgraph->template output<MeshData<MeshState::CorrectorPressure>>(parallelLane);
        subgraph->template input<MeshData<MeshState::CorrectorPressure>>(parallelLane);
    } else {
        auto corrPressureBarrier = makeTerminableRetaggingBarrier<
            MeshState::CorrectorPressure, MeshState::CorrectorPressure>(
            nmeshes, "CorrPressure",
            "PRESSURE_ITERATION",
            [](auto& meshes) {
                fds_pressure_iteration(meshes[0]->t, meshes[0]->dt);
            });
        subgraph->template edge<MeshData<MeshState::CorrectorPressure>>(
            parallelLane, corrPressureBarrier);
        subgraph->template edge<MeshData<MeshState::CorrectorPressure>>(
            corrPressureBarrier, parallelLane);
        subgraph->template input<TerminationData>(corrPressureBarrier);
    }

    // C7 → PostVelCorr → CorrFinalBarrier → PostCorrFinalBarrier → C8
    subgraph->template edge<MeshData<MeshState::PostVelCorr>>(
        parallelLane, corrFinalBarrier);
    subgraph->template edge<MeshData<MeshState::PostCorrFinalBarrier>>(
        corrFinalBarrier, parallelLane);

    // C8 → Default → subgraph output → TimestepTask
    subgraph->template output<MeshData<>>(parallelLane);

    // =====================================================================
    // TERMINATION WIRING — breaks all structural cycles
    // =====================================================================
    subgraph->template input<TerminationData>(forkableLane);
    subgraph->template input<TerminationData>(parallelLane);
    subgraph->template input<TerminationData>(hvacBarrier);
    subgraph->template input<TerminationData>(divExchangeTask);
    subgraph->template input<TerminationData>(fork2JoinTask);
    subgraph->template input<TerminationData>(corrFinalBarrier);

    return subgraph;
}

/// Dispatch wrapper: selects the correct template instantiation at runtime.
inline auto buildComputeSubgraph(int nmeshes, const ThreadBudget &budget) {
    if (fds_use_pressure_subgraph()) {
        return buildComputeSubgraphImpl<MeshState::PredictorPressure>(
            nmeshes, budget);
    }
    return buildComputeSubgraphImpl<MeshState::PredDivP2Out>(
        nmeshes, budget);
}

#endif // COMPUTE_SUBGRAPH_H
