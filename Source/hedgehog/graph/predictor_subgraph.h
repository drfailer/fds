#ifndef PREDICTOR_SUBGRAPH_H
#define PREDICTOR_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <service/comm_service.hpp>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/termination_data.h"
#include "../state/barrier_state.h"
#include "../task/pred_step1_kernel_task.h"
#include "../task/mass_fd_kernel_task.h"
#include "../task/div_setup_kernel_task.h"
#include "../task/pred_prefork_div_task.h"
#include "../task/divergence_part2_kernel_task.h"
#include "../task/velocity_predictor_kernel_task.h"
#include "../task/velocity_bc_edges_task.h"
#include "../task/barrier_tasks.h"
#include "../task/pred_div_parallel_task.h"
#include "../task/pred_cc_partmom_divp1_kernel_task.h"
#include "../tool/thread_budget.h"
#include "change_timestep_subgraph.h"
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
inline auto buildPredictorSubgraphImpl(int nmeshes, const ThreadBudget &budget,
                                    std::shared_ptr<MeshDependencyGraph> depGraph = nullptr,
                                    hh::comm::CommService *commService = nullptr) {
    auto subgraph = std::make_shared<hh::Graph<3,
        MeshData<>, TerminationData, MeshData<MeshState::PredictorPressure>,
        MeshData<>, MeshData<MeshState::PredictorPressure>>>("Predictor");

    // --- Kernel tasks (threads from budget) ---

    auto predStep1KernelTask = std::make_shared<PredStep1KernelTask>(budget.predStep1);
    auto velPredKernelTask = std::make_shared<VelocityPredictorKernelTask<PressureTag>>(budget.velPredictor);

    bool ccIBM = fds_is_cc_ibm() != 0;
    constexpr bool useParallelPressure = (PressureTag != MeshState::Default);

    auto changeTimeStepSubgraph = buildChangeTimeStepSubgraph(nmeshes, budget, ccIBM);

    // --- PredFinal: merged SynTurb+VelBC kernel → PhaseTransition ---
    auto predSynTurbVelBCTask = std::make_shared<PredSynTurbVelBCTask>(budget.predSynTurbVelBC);
    auto phaseTransTask = std::make_shared<PhaseTransitionTask>(nmeshes);

    // --- Common barrier states ---

    // (ChangeTimeStepCollector removed — RetryPreKernel collects N MeshData<> directly)

    // --- Wire the sub-graph ---

    subgraph->input<MeshData<>>(predStep1KernelTask);
    subgraph->input<TerminationData>(changeTimeStepSubgraph);

    // --- Predictor middle section: Fork (non-CC_IBM) or Sequential (CC_IBM) ---

    if (!ccIBM) {
        // Split barrier: MeshExchange(1) + Hvac + InitDiv (global only)
        // DivP1Prefork per-mesh loop extracted to downstream kernel task.
        auto meshExch1SM = makeBarrierSM(nmeshes, "MeshExch1+Hvac+InitDiv",
            "MESH_EXCHANGE(1)\\nEXCH_INS_PART\\nHVAC_CALC\\nINIT_DIV",
            [](auto& meshes) {
                fds_mesh_exchange(1);
                fds_exchange_inserted_particles();
                fds_hvac_calc(meshes[0]->t, meshes[0]->dt, 1);
                fds_initialize_divergence_integrals();
            });

        // Merged prefork + fork: DivP1Prefork + (DivSetup+PartMom || WallBC+DivP1Early)
        // AsyncWorker handles the fork internally — each Hedgehog thread owns a
        // worker thread, so real OS thread count = 2 × budget.predPreforkDiv.
        auto predPreforkDivTask = std::make_shared<PredPreforkDivTask>(budget.predPreforkDiv);

        subgraph->edges(predStep1KernelTask, meshExch1SM);
        subgraph->edges(meshExch1SM, predPreforkDivTask);

        // Divergence pipeline: packed parallel task + 2 retagging barriers.
        // The 3 parallel kernels (DivP1Late, DivP2Pre, DivPart2) share one
        // thread pool via PredDivParallelTask. Sequential barriers are separate
        // nodes, clearly highlighting the sequential/parallel structure.
        //
        // Pipeline: PredPreforkDiv → Parallel(DivP1Late) → DivExchange barrier
        //   → Parallel(DivP2Pre) → GlobalMatrix barrier → Parallel(DivPart2)
        //   → downstream

        auto predDivParallelTask = std::make_shared<PredDivParallelTask<PressureTag>>(
            budget.predDivParallel);

        auto divExchangeBarrier = makeRetaggingBarrier<MeshState::DivExch, MeshState::DivP2Pre>(
            nmeshes, "DivExchange",
            "EXCH_DIV_INFO",
            [](auto&) {
                fds_exchange_divergence_info();
            });

        auto globalMatBarrier = makeRetaggingBarrier<MeshState::GlobalMat, MeshState::DivPart2>(
            nmeshes, "GlobalMatrix+PressureInit",
            "GLOBAL_MATRIX_REASSIGN\\nPRESSURE_INIT",
            [useParallelPressure](auto&) {
                fds_global_matrix_reassign(0);
                if (useParallelPressure) {
                    fds_pressure_iteration_init();
                    fds_pressure_iteration_increment();
                }
            });

        // PredPreforkDiv → Parallel ↔ DivExchange ↔ Parallel ↔ GlobalMatrix ↔ Parallel
        subgraph->edges(predPreforkDivTask, predDivParallelTask);
        subgraph->edges(predDivParallelTask, divExchangeBarrier);
        subgraph->edges(divExchangeBarrier, predDivParallelTask);
        subgraph->edges(predDivParallelTask, globalMatBarrier);
        subgraph->edges(globalMatBarrier, predDivParallelTask);

        // TerminationData breaks structural cycle at shutdown
        subgraph->template input<TerminationData>(predDivParallelTask);

        // Downstream: final output (MeshData<PressureTag>) → Pressure → VelPred
        if constexpr (useParallelPressure) {
            subgraph->template output<MeshData<PressureTag>>(predDivParallelTask);
            subgraph->template input<MeshData<PressureTag>>(velPredKernelTask);
        } else {
            auto predPressureSM = makeBarrierSM(nmeshes, "PredPressure",
                "PRESSURE_ITERATION\\nINIT_CHANGE_TIME_STEP",
                [](auto& meshes) {
                    fds_pressure_iteration(meshes[0]->t, meshes[0]->dt);
                    fds_init_change_time_step(meshes[0]->dt);
                });
            subgraph->edges(predDivParallelTask, predPressureSM);
            subgraph->edges(predPressureSM, velPredKernelTask);
        }
    } else {
        // CC_IBM: MeshExchange(1) only
        auto meshExchange1SM = makeBarrierSM(nmeshes, "MeshExchange(1)",
            "MESH_EXCHANGE(1)\\nEXCHANGE_INSERTED_PARTICLES",
            [](auto& meshes) {
                fds_mesh_exchange(1);
                fds_exchange_inserted_particles();
            });

        subgraph->edges(predStep1KernelTask, meshExchange1SM);

        // CC_IBM sequential path
        auto hvacInitDivSM = makeBarrierSM(nmeshes, "Hvac+InitDiv",
            "HVAC_CALC\\nINITIALIZE_DIVERGENCE_INTEGRALS",
            [](auto& meshes) {
                fds_hvac_calc(meshes[0]->t, meshes[0]->dt, 1);
                fds_initialize_divergence_integrals();
            });

        auto predDivSetupKernelTask = std::make_shared<DivSetupKernelTask>(budget.standalone(4));
        subgraph->edges(meshExchange1SM, predDivSetupKernelTask);
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

    // VelocityPredictor -> ChangeTimeStep (RetryPreKernel collects N MeshData<> directly)
    subgraph->edges(velPredKernelTask, changeTimeStepSubgraph);

    // ChangeTimeStep scatters MeshData<> → SynTurb+VelBC(parallel) → PhaseTransition(collect+scatter)
    subgraph->edges(changeTimeStepSubgraph, predSynTurbVelBCTask);
    subgraph->edges(predSynTurbVelBCTask, phaseTransTask);

    subgraph->outputs(phaseTransTask);

    return subgraph;
}

/// Dispatch wrapper: selects the correct template instantiation at runtime.
inline auto buildPredictorSubgraph(int nmeshes, const ThreadBudget &budget,
                                    std::shared_ptr<MeshDependencyGraph> depGraph = nullptr,
                                    hh::comm::CommService *commService = nullptr) {
    if (fds_use_pressure_subgraph()) {
        return buildPredictorSubgraphImpl<MeshState::PredictorPressure>(
            nmeshes, budget, depGraph, commService);
    }
    return buildPredictorSubgraphImpl<MeshState::Default>(
        nmeshes, budget, depGraph, commService);
}

#endif // PREDICTOR_SUBGRAPH_H
