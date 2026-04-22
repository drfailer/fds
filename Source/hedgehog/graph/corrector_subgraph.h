#ifndef CORRECTOR_SUBGRAPH_H
#define CORRECTOR_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/termination_data.h"
#include "../state/barrier_state.h"
#include "../state/fork_join_state.h"
#include "../task/barrier_tasks.h"
#include "../task/mass_fd_kernel_task.h"
#include "../task/div_setup_kernel_task.h"
#include "../task/combustion_kernel_task.h"
#include "../task/corr_divsetup_comb_part_task.h"
#include "../task/corr_div_part1_kernel_task.h"
#include "../task/divergence_part2_kernel_task.h"
#include "../task/corr_final_kernel_task.h"
#include "../task/corr_div_parallel_task.h"
#include "../task/div_exchange_task.h"
#include "velocity_bc_subgraph.h"
#include "../task/corr_radiation_kernel_task.h"
#include "../task/pipeline_fork2_tasks.h"
#include "../tool/thread_budget.h"

/// Build the Corrector sub-graph.
///
/// Barrier splits applied:
///   - Group A: Split into pre-barrier (soot+hvac) + ParticleOpsKernel + post-barrier (exchange+WallBC)
///   - Group B: WallBCFinalize stays in barrier (uses POINT_TO_MESH, not thread-safe)
///   - Group C: Split into pre-barrier (MeshExch2) + QRAddCopyKernel + post-barrier (DivExch)
template<MeshState PressureTag = MeshState::Default>
inline auto buildCorrectorSubgraphImpl(int nmeshes, const ThreadBudget &budget) {
    auto subgraph = std::make_shared<hh::Graph<6,
        MeshData<>, MeshData<MeshState::CorrectorPressure>, MeshData<MeshState::PostCorrStep1>,
        MeshData<MeshState::PostParticleOps>, MeshData<MeshState::PostRadExch>, TerminationData,
        MeshData<>, MeshData<MeshState::CorrectorPressure>,
        MeshData<MeshState::MeshExch4>, MeshData<MeshState::MeshExch2>,
        MeshData<MeshState::MeshExch7>, BarrierData>>("Corrector");

    // --- Kernel tasks (threads from budget) ---

    auto corrFinalKernelTask = std::make_shared<CorrFinalKernelTask>(budget.corrFinal);

    bool ccIBM = fds_is_cc_ibm() != 0;
    bool ht3d = fds_is_ht3d() != 0;

    // --- Radiation kernel (replaces CorrRadiationSubgraph) ---
    auto corrRadiationKernelTask = std::make_shared<CorrRadiationKernelTask<MeshState::PostWallBC>>(
        budget.corrFork2Radiation);

    constexpr bool useParallelPressure = (PressureTag != MeshState::Default);

    // --- CorrStep1 + DivSetup||Comb + ParticleOps → ExchangeGraph → HvacCalc ---
    //
    // CorrDivSetupCombPartTask merges CorrStep1 + Fork1 (DivSetup||Comb via
    // AsyncWorker) + ParticleOps. Phase 1 (CorrStep1) → external exchange graph
    // → Phase 2:
    //   MeshData<> → HvacCalc barrier (emitted before ParticleOps)
    //   MeshData<MeshExch7> → exchange graph → PostParticleOps → HvacCalc
    // HvacCalc runs HVAC_CALC eagerly when N MeshData<> arrive, then waits
    // for N PostParticleOps before emitting N MeshData<PostHvac>.
    auto corrDivSetupCombPartTask = std::make_shared<CorrDivSetupCombPartTask>(
        budget.corrDivSetupCombPart);

    // HvacCalc: collects N MeshData<> + N MeshData<PostParticleOps>.
    // Eager: runs HVAC_CALC when primary (MeshData<>) set is complete.
    // Emits N MeshData<PostHvac> only when BOTH sets are complete.
    // Terminable: in cycle with CorrDivSetupCombPart (Phase 2 → HVAC → Phase 3).
    auto hvacBarrier = makeTerminableEagerDualMeshBarrier<
        MeshState::PostParticleOps, MeshState::PostHvac>(
        nmeshes, "HvacCalc",
        "HVAC_CALC",
        [](auto& meshes) {
            fds_hvac_calc(meshes[0]->t, meshes[0]->dt, 1);
        });

    // --- Wire the sub-graph ---

    // Subgraph input → CorrDivSetupCombPart (Phase 1: CorrStep1)
    subgraph->inputs(corrDivSetupCombPartTask);

    // CorrDivSetupCombPart Phase 1 → MeshExch4 → subgraph output (to exchange graph)
    // CC_IBM needs MESH_CC_EXCHANGE(4) before species exchange (cut-cell data sync)
    if (ccIBM) {
        auto ccExch4Barrier = makeBarrierSM<MeshState::MeshExch4>(
            nmeshes, "CC_Exchange(4)",
            "MESH_CC_EXCHANGE(4)",
            [](auto&) { fds_mesh_cc_exchange(4); });
        subgraph->edges(corrDivSetupCombPartTask, ccExch4Barrier);
        subgraph->template output<MeshData<MeshState::MeshExch4>>(ccExch4Barrier);
    } else {
        subgraph->template output<MeshData<MeshState::MeshExch4>>(corrDivSetupCombPartTask);
    }
    // Subgraph input (PostCorrStep1 from exchange graph) → CorrDivSetupCombPart Phase 2
    subgraph->template input<MeshData<MeshState::PostCorrStep1>>(corrDivSetupCombPartTask);

    // MeshData<> output → HvacCalc barrier (emitted before ParticleOps)
    subgraph->edges(corrDivSetupCombPartTask, hvacBarrier);
    // MeshData<MeshExch7> → subgraph output → exchange graph
    subgraph->template output<MeshData<MeshState::MeshExch7>>(corrDivSetupCombPartTask);
    // PostParticleOps from exchange graph → HvacCalc (secondary input)
    subgraph->template input<MeshData<MeshState::PostParticleOps>>(hvacBarrier);
    // HvacCalc (MeshData<PostHvac>) → CorrDivSetupCombPart Phase 3 (WallBC)
    subgraph->edges(hvacBarrier, corrDivSetupCombPartTask);

    // TerminationData breaks CorrDivSetupCombPart ↔ HvacCalc cycle
    subgraph->template input<TerminationData>(hvacBarrier);

    // --- Fork 2: RADIATION || DIV_P1 (or sequential for CC_IBM) ---

    // Join barrier: collects N MeshData<> (from DivP1 path) + N MeshData<PostRadExch>
    // (from radiation exchange graph), emits N MeshData<>.
    auto fork2JoinTask = std::make_shared<DualMeshJoinTask<MeshState::PostRadExch>>(
        nmeshes, "Fork2Join");

    // Radiation → MeshExch2 → exchange graph → PostRadExch → join
    subgraph->edges(corrDivSetupCombPartTask, corrRadiationKernelTask);
    subgraph->template output<MeshData<MeshState::MeshExch2>>(corrRadiationKernelTask);
    subgraph->template input<MeshData<MeshState::PostRadExch>>(fork2JoinTask);

    if (ccIBM) {
        // Group B (CC_IBM): Exchange(6) barrier for back wall data (HT3D only) + InitDiv.
        auto groupBPostSM = makeRetaggingBarrier<MeshState::PostWallBC, MeshState::Default>(
            nmeshes, "MeshExch6a+InitDiv",
            "MESH_EXCHANGE(6) [HT3D]\\nINIT_DIV_INTEGRALS",
            [ht3d](auto& meshes) {
                if (ht3d && meshes[0]->call_ht_1d) { fds_mesh_exchange(6); }
                fds_initialize_divergence_integrals();
            });

        auto corrDivP1KernelTask = std::make_shared<CorrDivPart1KernelTask>(budget.standalone(2));

        // CC_IBM: per-mesh loop kept for GET_LINKED_VELOCITIES (cross-mesh writes)
        auto corrDivExchangeSM = makeBarrierSM(nmeshes, "CorrDivExchange",
            "EXCH_DIV_INFO\\nZONE_OPS\\nGLOBAL_MATRIX_REASSIGN\\nPRES_INIT+INCR",
            [useParallelPressure](auto& meshes) {
                fds_exchange_divergence_info();
                for (auto &md : meshes) {
                    fds_divergence_part_2_preprocessing(md->nm, md->dt);
                }
                fds_global_matrix_reassign(0);
                if (useParallelPressure) {
                    fds_pressure_iteration_init();
                    fds_pressure_iteration_increment();
                    for (auto &md : meshes) {
                        fds_get_linked_fv(md->nm, 0);
                    }
                }
            });

        // Fork: Radiation starts immediately, Exchange(6)+InitDiv runs in parallel.
        // Both join at fork2JoinTask.
        subgraph->edges(corrDivSetupCombPartTask, groupBPostSM);
        subgraph->edges(groupBPostSM, fork2JoinTask);
        subgraph->edges(fork2JoinTask, corrDivP1KernelTask);
        subgraph->edges(corrDivP1KernelTask, corrDivExchangeSM);
        // CC_IBM keeps standalone DivPart2 kernel (not packed)
        auto corrDivP2KernelTask = std::make_shared<DivergencePart2KernelTask<PressureTag>>(budget.corrDivPart2);
        subgraph->edges(corrDivExchangeSM, corrDivP2KernelTask);

        // Downstream: DivP2 → Pressure → CorrFinalKernel
        if constexpr (useParallelPressure) {
            subgraph->outputs(corrDivP2KernelTask);
            subgraph->template input<MeshData<MeshState::CorrectorPressure>>(corrFinalKernelTask);
        } else {
            auto corrPressureSM = makeRetaggingBarrier<MeshState::Default, MeshState::CorrectorPressure>(
                nmeshes, "CorrPressure",
                "PRESSURE_ITERATION",
                [](auto& meshes) {
                    fds_pressure_iteration(meshes[0]->t, meshes[0]->dt);
                });
            subgraph->edges(corrDivP2KernelTask, corrPressureSM);
            subgraph->edges(corrPressureSM, corrFinalKernelTask);
        }
    } else {
        // Group B (non-CC_IBM): Exchange(6) [HT3D only] + InitDiv barrier.
        auto groupBPostSM = makeRetaggingBarrier<MeshState::PostWallBC, MeshState::Default>(
            nmeshes, "MeshExch6a+InitDiv",
            "MESH_EXCHANGE(6) [HT3D]\\nINIT_DIV_INTEGRALS",
            [ht3d](auto& meshes) {
                if (ht3d && meshes[0]->call_ht_1d) { fds_mesh_exchange(6); }
                fds_initialize_divergence_integrals();
            });

        auto fork2DivP1Task = std::make_shared<Fork2DivP1KernelTask>(budget.corrFork2DivP1);

        // Divergence pipeline: packed parallel task + DivExchangeTask.
        auto corrDivParallelTask = std::make_shared<CorrDivParallelTask<PressureTag>>(
            budget.corrDivParallel);

        auto corrDivExchangeTask = std::make_shared<DivExchangeTask<PressureTag>>(
            nmeshes, budget.divExchange);

        // Fork: Radiation starts immediately.
        //        DivP1 waits for Exchange(6)+InitDiv barrier.
        subgraph->edges(corrDivSetupCombPartTask, groupBPostSM);
        subgraph->edges(groupBPostSM, fork2DivP1Task);

        // Both branches join at fork2JoinTask.
        subgraph->edges(fork2DivP1Task, fork2JoinTask);
        subgraph->edges(fork2JoinTask, corrDivParallelTask);

        // CorrDivParallel ↔ DivExchangeTask
        subgraph->edges(corrDivParallelTask, corrDivExchangeTask);
        subgraph->edges(corrDivExchangeTask, corrDivParallelTask);

        // TerminationData breaks structural cycle at shutdown
        subgraph->template input<TerminationData>(corrDivExchangeTask);

        // Downstream: final output (MeshData<PressureTag>) → Pressure → CorrFinalKernel
        if constexpr (useParallelPressure) {
            subgraph->template output<MeshData<PressureTag>>(corrDivParallelTask);
            subgraph->template input<MeshData<MeshState::CorrectorPressure>>(corrFinalKernelTask);
        } else {
            auto corrPressureSM = makeRetaggingBarrier<MeshState::Default, MeshState::CorrectorPressure>(
                nmeshes, "CorrPressure",
                "PRESSURE_ITERATION",
                [](auto& meshes) {
                    fds_pressure_iteration(meshes[0]->t, meshes[0]->dt);
                });
            subgraph->edges(corrDivParallelTask, corrPressureSM);
            subgraph->edges(corrPressureSM, corrFinalKernelTask);
        }
    }

    // --- CorrFinal: CorrFinalKernel ↔ CorrFinalOrch → CorrFinalDump ---
    //
    // CorrFinalKernelTask has 3 phases via different input types:
    //   Phase 1: MeshData<CorrectorPressure> → VelCorr → MeshData<PostVelCorr> → Orch
    //   Phase 2: MeshData<> (from Orch) → VelBCEdges → MeshData<> → Dump
    //   Phase 3: BarrierData (from Orch) → RTE → BarrierData → Dump

    auto corrFinalOrchTask = std::make_shared<CorrFinalOrchTask>(nmeshes, ccIBM);
    auto corrFinalDumpTask = std::make_shared<CorrFinalDumpTask>(nmeshes);

    // CorrFinalKernel ↔ CorrFinalOrch cycle
    subgraph->edges(corrFinalKernelTask, corrFinalOrchTask);  // MeshData<PostVelCorr>
    subgraph->edges(corrFinalOrchTask, corrFinalKernelTask);  // MeshData<> + BarrierData
    // CorrFinalKernel → CorrFinalDump
    subgraph->edges(corrFinalKernelTask, corrFinalDumpTask);  // MeshData<> + BarrierData

    // TerminationData breaks CorrFinalKernel ↔ CorrFinalOrch cycle
    subgraph->template input<TerminationData>(corrFinalOrchTask);

    subgraph->outputs(corrFinalDumpTask);

    return subgraph;
}

/// Dispatch wrapper: selects the correct template instantiation at runtime.
inline auto buildCorrectorSubgraph(int nmeshes, const ThreadBudget &budget) {
    if (fds_use_pressure_subgraph()) {
        return buildCorrectorSubgraphImpl<MeshState::CorrectorPressure>(
            nmeshes, budget);
    }
    return buildCorrectorSubgraphImpl<MeshState::Default>(
        nmeshes, budget);
}

#endif // CORRECTOR_SUBGRAPH_H
