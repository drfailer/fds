#ifndef PREDICTOR_SUBGRAPH_H
#define PREDICTOR_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/termination_data.h"
#include "../state/barrier_state.h"
#include "../task/pred_step1_kernel_task.h"
#include "../task/mass_fd_kernel_task.h"
#include "../task/div_setup_kernel_task.h"
#include "../task/pred_prefork_div_task.h"
#include "../task/div_exchange_task.h"
#include "../task/divergence_part2_kernel_task.h"
#include "../task/velocity_predictor_kernel_task.h"
#include "../task/velocity_bc_edges_task.h"
#include "../task/barrier_tasks.h"
#include "../task/pred_cc_partmom_divp1_kernel_task.h"
#include "../task/change_timestep_task.h"
#include "../tool/thread_budget.h"
#include "velocity_bc_subgraph.h"
#include "../task/wallbc_kernel_task.h"

/// Build the Predictor sub-graph.
///
/// Template parameter PressureTag controls the DivP2 output and VelPred input
/// MeshData state. When PressureTag != Default, DivP2 outputs typed data that
/// exits the subgraph for the shared pressure subgraph, and VelPred accepts
/// the typed result back.
///
/// Barrier splits applied:
///   - meshExch1DivPrefork: DivP1Prefork per-mesh loop → parallel kernel task
///   - meshExch3SynTurb: SyntheticTurbulence per-mesh loop → parallel kernel task
///   - predJoinDivExchange: DivP1Late per-mesh loop → ForkJoin + parallel kernel task
template<MeshState PressureTag = MeshState::Default>
inline auto buildPredictorSubgraphImpl(int nmeshes, const ThreadBudget &budget) {
    auto subgraph = std::make_shared<hh::Graph<5,
        MeshData<>, TerminationData, MeshData<MeshState::PredictorPressure>,
        MeshData<MeshState::PostPredExch>, MeshData<MeshState::PostPredVelExch>,
        MeshData<>, MeshData<MeshState::PredictorPressure>,
        MeshData<MeshState::MeshExch1>, MeshData<MeshState::MeshExch3>>>("Predictor");

    // --- Kernel tasks (threads from budget) ---

    auto predStep1KernelTask = std::make_shared<PredStep1KernelTask>(budget.predStep1);
    auto velPredKernelTask = std::make_shared<VelocityPredictorKernelTask<PressureTag>>(budget.velPredictor);

    bool ccIBM = fds_is_cc_ibm() != 0;
    constexpr bool useParallelPressure = (PressureTag != MeshState::Default);

    auto changeTimeStepTask = std::make_shared<ChangeTimeStepTask>(
        nmeshes, budget.retryMomDiv, ccIBM);

    // --- PredFinal: merged SynTurb+VelBC kernel → PhaseTransition ---
    auto predSynTurbVelBCTask = std::make_shared<PredSynTurbVelBCTask>(budget.predSynTurbVelBC);
    auto phaseTransTask = std::make_shared<PhaseTransitionTask>(nmeshes);

    // --- Common barrier states ---

    // --- Wire the sub-graph ---

    subgraph->input<MeshData<>>(predStep1KernelTask);
    // ChangeTimeStepTask is a single task (no cycle) — no TerminationData needed

    // --- Predictor middle section: Fork (non-CC_IBM) or Sequential (CC_IBM) ---

    if (!ccIBM) {
        // MeshExchange(1) routed through exchange graph.
        // Post-exchange barrier handles remaining global ops.
        auto postPredExchBarrier = makeRetaggingBarrier<MeshState::PostPredExch, MeshState::Default>(
            nmeshes, "ExchInsPart+Hvac+InitDiv",
            "EXCH_INS_PART\\nHVAC_CALC\\nINIT_DIV",
            [](auto& meshes) {
                fds_exchange_inserted_particles();
                fds_hvac_calc(meshes[0]->t, meshes[0]->dt, 1);
                fds_initialize_divergence_integrals();
            });

        // Merged prefork + fork + divergence pipeline: DivP1Prefork +
        // (DivSetup+PartMom || WallBC+DivP1Early) + DivP1Late + DivPart2.
        // AsyncWorker handles the fork internally — each Hedgehog thread owns a
        // worker thread, so real OS thread count = 2 × budget.predPreforkDiv.
        //
        // Pipeline: PredPreforkDiv(Phase1) → DivExchangeTask (collect+parallel DivP2Pre)
        //   → PredPreforkDiv(Phase2) → downstream
        auto predPreforkDivTask = std::make_shared<PredPreforkDivTask<PressureTag>>(
            budget.predPreforkDiv);

        auto divExchangeTask = std::make_shared<DivExchangeTask<PressureTag>>(
            nmeshes, budget.divExchange);

        subgraph->template output<MeshData<MeshState::MeshExch1>>(predStep1KernelTask);
        subgraph->template input<MeshData<MeshState::PostPredExch>>(postPredExchBarrier);
        subgraph->edges(postPredExchBarrier, predPreforkDivTask);

        // PredPreforkDiv ↔ DivExchangeTask
        subgraph->edges(predPreforkDivTask, divExchangeTask);
        subgraph->edges(divExchangeTask, predPreforkDivTask);

        // TerminationData breaks structural cycle at shutdown
        subgraph->template input<TerminationData>(divExchangeTask);

        // Downstream: final output (MeshData<PressureTag>) → Pressure → VelPred
        if constexpr (useParallelPressure) {
            subgraph->template output<MeshData<PressureTag>>(predPreforkDivTask);
            subgraph->template input<MeshData<PressureTag>>(velPredKernelTask);
        } else {
            auto predPressureSM = makeBarrierSM(nmeshes, "PredPressure",
                "PRESSURE_ITERATION\\nINIT_CHANGE_TIME_STEP",
                [](auto& meshes) {
                    fds_pressure_iteration(meshes[0]->t, meshes[0]->dt);
                    fds_init_change_time_step(meshes[0]->dt);
                });
            subgraph->edges(predPreforkDivTask, predPressureSM);
            subgraph->edges(predPressureSM, velPredKernelTask);
        }
    } else {
        // CC_IBM: MeshExchange(1) routed through exchange graph.
        // Post-exchange barrier handles exchange_inserted_particles only.
        auto postPredExchBarrier = makeRetaggingBarrier<MeshState::PostPredExch, MeshState::Default>(
            nmeshes, "ExchInsPart",
            "EXCH_INS_PART",
            [](auto&) {
                fds_exchange_inserted_particles();
            });

        subgraph->template output<MeshData<MeshState::MeshExch1>>(predStep1KernelTask);
        subgraph->template input<MeshData<MeshState::PostPredExch>>(postPredExchBarrier);

        // CC_IBM sequential path
        auto hvacInitDivSM = makeBarrierSM(nmeshes, "Hvac+InitDiv",
            "HVAC_CALC\\nINITIALIZE_DIVERGENCE_INTEGRALS",
            [](auto& meshes) {
                fds_hvac_calc(meshes[0]->t, meshes[0]->dt, 1);
                fds_initialize_divergence_integrals();
            });

        auto predDivSetupKernelTask = std::make_shared<DivSetupKernelTask>(budget.standalone(4));
        subgraph->edges(postPredExchBarrier, predDivSetupKernelTask);
        subgraph->edges(predDivSetupKernelTask, hvacInitDivSM);

        // WallBC inlined: no orchestrator needed in predictor (dt_bc=0, call_ht_1d=0 defaults)
        // WallBCKernelTask includes finalize (all per-mesh, thread-safe)
        auto predWallBCKernel = std::make_shared<WallBCKernelTask>(budget.standalone(2));

        // Extracted: PartMom + DivP1 per mesh (parallel, Loop 1 from original barrier)
        auto predCCPartMomDivP1Kernel = std::make_shared<PredCCPartMomDivP1KernelTask>(budget.standalone(2));

        // Barrier shrunk: Loop 1 extracted, only exchange + Loop 2 + global ops remain
        // CC_IBM: Loop 2 per-mesh loop kept for GET_LINKED_VELOCITIES (cross-mesh writes)
        auto predDivExchangeSM = makeBarrierSM(nmeshes, "DivExch+ZoneOps",
            "EXCH_DIV_INFO\\nZONE_OPS\\nGLOBAL_MATRIX_REASSIGN\\nPRES_INIT+INCR",
            [useParallelPressure](auto& meshes) {
                fds_exchange_divergence_info();
                // Zone ops + GET_LINKED_VELOCITIES (CC_IBM needs per-mesh for cross-mesh writes)
                for (auto &md : meshes) {
                    fds_divergence_part_2_preprocessing(md->nm, md->dt);
                }
                fds_global_matrix_reassign(0);
                if (useParallelPressure) {
                    fds_pressure_iteration_init();
                    fds_pressure_iteration_increment();
                    // Pre-loop: link cut-face velocity fluxes before pressure iterations
                    for (auto &md : meshes) {
                        fds_get_linked_fv(md->nm, 0); // DO_BAROCLINIC=FALSE
                    }
                }
            });

        subgraph->edges(hvacInitDivSM, predWallBCKernel);
        subgraph->edges(predWallBCKernel, predCCPartMomDivP1Kernel);
        subgraph->edges(predCCPartMomDivP1Kernel, predDivExchangeSM);

        // CC_IBM keeps standalone DivPart2 kernel (not packed)
        auto predDivP2KernelTask = std::make_shared<DivergencePart2KernelTask<PressureTag>>(budget.predDivPart2);
        subgraph->edges(predDivExchangeSM, predDivP2KernelTask);

        // Downstream: DivP2 → Pressure → VelPred
        if constexpr (useParallelPressure) {
            subgraph->outputs(predDivP2KernelTask);
            subgraph->template input<MeshData<PressureTag>>(velPredKernelTask);
        } else {
            auto predPressureSM = makeBarrierSM(nmeshes, "PredPressure",
                "PRESSURE_ITERATION\\nINIT_CHANGE_TIME_STEP",
                [](auto& meshes) {
                    fds_pressure_iteration(meshes[0]->t, meshes[0]->dt);
                    fds_init_change_time_step(meshes[0]->dt);
                });
            subgraph->edges(predDivP2KernelTask, predPressureSM);
            subgraph->edges(predPressureSM, velPredKernelTask);
        }
    }

    // VelocityPredictor → ChangeTimeStep (collects N, retries internally)
    subgraph->edges(velPredKernelTask, changeTimeStepTask);

    // ChangeTimeStep → MeshExch3 → exchange graph → PostPredVelExch → SynTurb+VelBC
    subgraph->template output<MeshData<MeshState::MeshExch3>>(changeTimeStepTask);
    subgraph->template input<MeshData<MeshState::PostPredVelExch>>(predSynTurbVelBCTask);
    subgraph->edges(predSynTurbVelBCTask, phaseTransTask);

    subgraph->outputs(phaseTransTask);

    return subgraph;
}

/// Dispatch wrapper: selects the correct template instantiation at runtime.
inline auto buildPredictorSubgraph(int nmeshes, const ThreadBudget &budget) {
    if (fds_use_pressure_subgraph()) {
        return buildPredictorSubgraphImpl<MeshState::PredictorPressure>(
            nmeshes, budget);
    }
    return buildPredictorSubgraphImpl<MeshState::Default>(
        nmeshes, budget);
}

#endif // PREDICTOR_SUBGRAPH_H
