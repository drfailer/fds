#ifndef PREDICTOR_SUBGRAPH_H
#define PREDICTOR_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <service/comm_service.hpp>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/termination_data.h"
#include "../state/collector_state.h"
#include "../state/barrier_state.h"
#include "../state/pred_step1_state.h"
#include "../state/div_setup_state.h"
#include "../state/fork_join_state.h"
#include "../task/pred_step1_kernel_task.h"
#include "../task/mass_fd_kernel_task.h"
#include "../task/div_setup_kernel_task.h"
#include "../task/pred_fork_tasks.h"
#include "../task/pred_wall_div_kernel_task.h"
#include "../task/divergence_part2_kernel_task.h"
#include "../task/velocity_predictor_kernel_task.h"
#include "../task/div_p1_prefork_kernel_task.h"
#include "../task/synthetic_turbulence_kernel_task.h"
#include "../task/div_p1_late_kernel_task.h"
#include "../tool/thread_budget.h"
#include "change_timestep_subgraph.h"
#include "velocity_bc_subgraph.h"
#include "../task/wallbc_kernel_task.h"
#include "pred_fork_vflux_subgraph.h"
#include "pred_fork_div_subgraph.h"

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

    auto predStep1OrchTask = std::make_shared<PredStep1Orchestrator>(nmeshes);
    auto predStep1KernelTask = std::make_shared<PredStep1KernelTask>(budget.predStep1);
    auto predDivP2KernelTask = std::make_shared<DivergencePart2KernelTask<PressureTag>>(budget.predDivPart2);
    auto velPredKernelTask = std::make_shared<VelocityPredictorKernelTask<PressureTag>>(budget.velPredictor);

    bool ccIBM = fds_is_cc_ibm() != 0;
    constexpr bool useParallelPressure = (PressureTag != MeshState::Default);

    auto predFinalSubgraph = buildPredFinalSubgraph(nmeshes, budget.predFinalVelBC);
    auto changeTimeStepSubgraph = buildChangeTimeStepSubgraph(nmeshes, budget.retryMomDiv);

    // --- Common barrier states ---

    auto changeTimeStepCollectorTask = std::make_shared<CollectorTask>(nmeshes, "ChangeTimeStepCollector");

    // Split barrier: MeshExchange(3) + CC_END_STEP (global only)
    // SyntheticTurbulence per-mesh loop extracted to downstream kernel task.
    auto meshExch3Task = makeBarrierTask("MeshExch3",
        "CC_END_STEP\\nMESH_EXCHANGE(3)",
        [ccIBM](auto& meshes) {
            if (ccIBM) { fds_cc_end_step(meshes[0]->t, meshes[0]->dt, 0); }
            fds_mesh_exchange(3);
        });

    // Extracted: SyntheticTurbulence per-mesh (parallel)
    auto synTurbKernelTask = std::make_shared<SyntheticTurbulenceKernelTask>(budget.predSynTurb);

    // --- Wire the sub-graph ---

    subgraph->input<MeshData<>>(predStep1OrchTask);
    subgraph->input<TerminationData>(changeTimeStepSubgraph);

    // PredStep1 (VISC + MASS_FD + DENSITY merged)
    subgraph->edges(predStep1OrchTask, predStep1KernelTask);

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

        // Extracted: DivP1Prefork per-mesh (parallel)
        auto divP1PreforkKernelTask = std::make_shared<DivP1PreforkKernelTask>(budget.predDivPrefork);

        subgraph->edges(predStep1KernelTask, meshExch1SM);
        subgraph->edges(meshExch1SM, divP1PreforkKernelTask);

        // Fork: (VFLUX + PART_MOM) || (WallBC + DIV_P1_early) — threads from budget
        auto predForkVFluxSG = buildPredForkVFluxSubgraph(
            nmeshes, budget.predForkDivSetup, budget.predForkPartMom);
        auto predForkDivSG = buildPredForkDivSubgraph(
            nmeshes, budget.predForkWallBC, budget.predForkDivP1Early);

        subgraph->edges(divP1PreforkKernelTask, predForkVFluxSG);
        subgraph->edges(divP1PreforkKernelTask, predForkDivSG);

        // Split barrier: ForkJoin(2N→N) → DivP1Late(parallel) → DivExchange(global)
        auto predForkJoinTask = std::make_shared<ForkJoinTask>(nmeshes, 2, budget.predDivP1Late, "PredForkJoin");

        // Extracted: DivP1Late per-mesh (parallel)
        auto divP1LateKernelTask = std::make_shared<DivP1LateKernelTask>(budget.predDivP1Late);

        // Smaller barrier: global divergence exchange + zone ops (no per-mesh DivP1Late)
        auto predDivExchangeSM = makeBarrierSM(nmeshes, "DivExch+DivP2Preproc",
            "EXCH_DIV_INFO\\nDIV_P2_PREPROC\\nGLOBAL_MATRIX_REASSIGN\\nPRES_INIT+INCR",
            [useParallelPressure](auto& meshes) {
                fds_exchange_divergence_info();
                // Zone ops for DivP2 (modifies global USUM, must run sequentially)
                for (auto &md : meshes) {
                    fds_divergence_part_2_preprocessing(md->nm, md->dt);
                }
                fds_global_matrix_reassign(0);
                if (useParallelPressure) {
                    fds_pressure_iteration_init();
                    fds_pressure_iteration_increment();
                }
            });

        subgraph->edges(predForkVFluxSG, predForkJoinTask);
        subgraph->edges(predForkDivSG, predForkJoinTask);
        subgraph->edges(predForkJoinTask, divP1LateKernelTask);
        subgraph->edges(divP1LateKernelTask, predDivExchangeSM);

        // DivP2 -> Pressure -> VelPred
        subgraph->edges(predDivExchangeSM, predDivP2KernelTask);
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

        auto predDivSetupOrchTask = std::make_shared<PredDivSetupOrchestrator>(nmeshes);
        auto predDivSetupKernelTask = std::make_shared<DivSetupKernelTask>(budget.standalone(4));
        subgraph->edges(meshExchange1SM, predDivSetupOrchTask);
        subgraph->edges(predDivSetupOrchTask, predDivSetupKernelTask);
        subgraph->edges(predDivSetupKernelTask, hvacInitDivSM);

        // WallBC inlined: no orchestrator needed in predictor (dt_bc=0, call_ht_1d=0 defaults)
        // WallBCKernelTask includes finalize (all per-mesh, thread-safe)
        auto predWallBCKernel = std::make_shared<WallBCKernelTask>(budget.standalone(2));

        // Barrier: PartMom + DivP1 + DivExchange + PressureInit
        auto predDivExchangeSM = makeBarrierSM(nmeshes, "WallDiv+DivExch",
            "PART_MOM\\nDIV_P1\\nEXCH_DIV_INFO\\nDIV_P2_PREPROC\\nGLOBAL_MATRIX_REASSIGN\\nPRES_INIT+INCR",
            [useParallelPressure](auto& meshes) {
                for (auto &md : meshes) {
                    fds_particle_momentum_kernel(md->nm, md->dt);
                    fds_divergence_part_1_kernel(md->nm, md->t, md->dt);
                }
                fds_exchange_divergence_info();
                // Zone ops for DivP2 (modifies global USUM, must run sequentially)
                for (auto &md : meshes) {
                    fds_divergence_part_2_preprocessing(md->nm, md->dt);
                }
                fds_global_matrix_reassign(0);
                if (useParallelPressure) {
                    fds_pressure_iteration_init();
                    fds_pressure_iteration_increment();
                }
            });

        subgraph->edges(hvacInitDivSM, predWallBCKernel);
        subgraph->edges(predWallBCKernel, predDivExchangeSM);

        // DivP2 -> Pressure -> VelPred
        subgraph->edges(predDivExchangeSM, predDivP2KernelTask);
    }

    // --- Common downstream: DivP2 -> Pressure -> VelPred -> ... ---

    if constexpr (useParallelPressure) {
        // Pressure handled externally via shared pressure subgraph.
        // DivP2<PressureTag> exits subgraph, VelPred<PressureTag> receives from outside.
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

    // VelocityPredictor -> ChangeTimeStep
    subgraph->edges(velPredKernelTask, changeTimeStepCollectorTask);
    subgraph->edges(changeTimeStepCollectorTask, changeTimeStepSubgraph);

    // Split: MeshExch(3) → SyntheticTurbulence(parallel) → PredFinal
    subgraph->edges(changeTimeStepSubgraph, meshExch3Task);
    subgraph->edges(meshExch3Task, synTurbKernelTask);
    subgraph->edges(synTurbKernelTask, predFinalSubgraph);

    subgraph->outputs(predFinalSubgraph);

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
