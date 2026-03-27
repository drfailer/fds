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
#include "../task/pred_step1_kernel_task.h"
#include "../task/mass_fd_kernel_task.h"
#include "../task/div_setup_kernel_task.h"
#include "../task/pred_fork_tasks.h"
#include "../task/pred_wall_div_kernel_task.h"
#include "../task/divergence_part2_kernel_task.h"
#include "../task/velocity_predictor_kernel_task.h"
#include "../tool/thread_budget.h"
#include "change_timestep_subgraph.h"
#include "velocity_bc_subgraph.h"
#include "../task/wallbc_kernel_task.h"
#include "pressure_iteration_subgraph.h"
#include "pred_fork_vflux_subgraph.h"
#include "pred_fork_div_subgraph.h"

/// Build the Predictor sub-graph.
///
/// Optimizations applied:
///   - MeshExchange(3) + PredFinalOrch merged into single barrier
///   - PredFinalCollector + PhaseTransition merged into PredFinal collector
///   - PredFinal subgraph now outputs MeshData directly (no BarrierData)
inline auto buildPredictorSubgraph(int nmeshes, const ThreadBudget &budget,
                                    std::shared_ptr<MeshDependencyGraph> depGraph = nullptr,
                                    hh::comm::CommService *commService = nullptr) {
    auto subgraph = std::make_shared<hh::Graph<2, MeshData, TerminationData, MeshData>>("Predictor");

    // --- Kernel tasks (threads from budget) ---

    auto predStep1OrchSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<PredStep1Orchestrator>(nmeshes), "PredStep1Orch");
    auto predStep1KernelTask = std::make_shared<PredStep1KernelTask>(budget.predStep1);
    auto predDivP2KernelTask = std::make_shared<DivergencePart2KernelTask>(budget.predDivPart2);
    auto velPredKernelTask = std::make_shared<VelocityPredictorKernelTask>(budget.velPredictor);

    bool ccIBM = fds_is_cc_ibm() != 0;
    bool useParallelPressure = fds_use_pressure_subgraph() != 0;

    auto predFinalSubgraph = buildPredFinalSubgraph(nmeshes, budget.predFinalVelBC);
    auto changeTimeStepSubgraph = buildChangeTimeStepSubgraph(nmeshes, budget.retryMomDiv);

    // --- Common barrier states ---

    auto changeTimeStepCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "ChangeTimeStepCollector");

    // Merged: MeshExchange(3) + PredFinalOrch (synthetic turbulence)
    auto meshExch3SynTurbSM = makeBarrierSM(nmeshes, "MeshExch3+SynTurb",
        "CC_END_STEP\\nMESH_EXCHANGE(3)\\nSYNTHETIC_TURBULENCE",
        [ccIBM](auto& meshes) {
            if (ccIBM) { fds_cc_end_step(meshes[0]->t, meshes[0]->dt, 0); }
            fds_mesh_exchange(3);
            for (auto &md : meshes) {
                fds_synthetic_turbulence_if_enabled(md->dt, md->t, md->nm);
            }
        });

    // --- Wire the sub-graph ---

    subgraph->input<MeshData>(predStep1OrchSM);
    subgraph->input<TerminationData>(changeTimeStepSubgraph);

    // PredStep1 (VISC + MASS_FD + DENSITY merged)
    subgraph->edges(predStep1OrchSM, predStep1KernelTask);

    // --- Predictor middle section: Fork (non-CC_IBM) or Sequential (CC_IBM) ---

    if (!ccIBM) {
        // Merged barrier: MeshExchange(1) + HvacInitDivPrefork
        auto meshExch1DivPreforkSM = makeBarrierSM(nmeshes, "MeshExch1+DivPrefork",
            "MESH_EXCHANGE(1)\\nEXCH_INS_PART\\nHVAC_CALC\\nINIT_DIV\\nDIV_P1_PREFORK",
            [](auto& meshes) {
                fds_mesh_exchange(1);
                fds_exchange_inserted_particles();
                fds_hvac_calc(meshes[0]->t, meshes[0]->dt, 1);
                fds_initialize_divergence_integrals();
                for (auto &md : meshes) {
                    fds_divergence_part_1_prefork(md->nm, md->t, md->dt);
                }
            });

        subgraph->edges(predStep1KernelTask, meshExch1DivPreforkSM);

        // Fork: (VFLUX + PART_MOM) || (WallBC + DIV_P1_early) — threads from budget
        auto predForkVFluxSG = buildPredForkVFluxSubgraph(
            nmeshes, budget.predForkDivSetup, budget.predForkPartMom);
        auto predForkDivSG = buildPredForkDivSubgraph(
            nmeshes, budget.predForkWallBC, budget.predForkDivP1Early);

        subgraph->edges(meshExch1DivPreforkSM, predForkVFluxSG);
        subgraph->edges(meshExch1DivPreforkSM, predForkDivSG);

        // Merged barrier: PredJoin + DivP1Late + PredDivExchange + PressureInit
        auto predJoinDivExchangeSM = makeBarrierSM(nmeshes, "PredJoin+DivLate+DivExch",
            "DIV_P1_LATE\\nEXCH_DIV_INFO\\nGLOBAL_MATRIX_REASSIGN\\nPRES_INIT+INCR",
            [useParallelPressure](auto& meshes) {
                for (auto &md : meshes) {
                    fds_divergence_part_1_late_b(md->nm, md->t, md->dt);
                }
                fds_exchange_divergence_info();
                fds_global_matrix_reassign(0);
                if (useParallelPressure) {
                    fds_pressure_iteration_init();
                    fds_pressure_iteration_increment();
                }
            },
            2 * nmeshes);  // Expects 2*N tokens from 2 fork branches

        subgraph->edges(predForkVFluxSG, predJoinDivExchangeSM);
        subgraph->edges(predForkDivSG, predJoinDivExchangeSM);

        // DivP2 -> Pressure -> VelPred
        subgraph->edges(predJoinDivExchangeSM, predDivP2KernelTask);
    } else {
        // CC_IBM: MeshExchange(1) only
        auto meshExchange1SM = makeBarrierSM(nmeshes, "MeshExchange(1)",
            "CC_DENSITY\\nMESH_EXCHANGE(1)\\nEXCHANGE_INSERTED_PARTICLES",
            [](auto& meshes) {
                fds_cc_density(meshes[0]->t, meshes[0]->dt);
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

        auto predDivSetupOrchSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
            std::make_shared<PredDivSetupOrchestrator>(nmeshes), "PredDivSetupOrch");
        auto predDivSetupKernelTask = std::make_shared<DivSetupKernelTask>(budget.standalone(4));
        subgraph->edges(meshExchange1SM, predDivSetupOrchSM);
        subgraph->edges(predDivSetupOrchSM, predDivSetupKernelTask);
        subgraph->edges(predDivSetupKernelTask, hvacInitDivSM);

        // WallBC inlined: no orchestrator needed in predictor (dt_bc=0, call_ht_1d=0 defaults)
        auto predWallBCKernel = std::make_shared<WallBCKernelTask>(budget.standalone(2));

        // Merged barrier: WallBCFinalize + PredWallDivKernel + PredDivExchange + PressureInit
        auto predDivExchangeSM = makeBarrierSM(nmeshes, "WallBCFin+WallDiv+DivExch",
            "WALLBC_FINALIZE\\nPART_MOM\\nDIV_P1\\nEXCH_DIV_INFO\\nGLOBAL_MATRIX_REASSIGN\\nPRES_INIT+INCR",
            [useParallelPressure](auto& meshes) {
                for (auto &md : meshes) {
                    fds_wall_bc_finalize(md->nm, md->t, md->dt_bc, md->call_ht_1d);
                }
                for (auto &md : meshes) {
                    fds_particle_momentum_kernel(md->nm, md->dt);
                    fds_divergence_part_1_kernel(md->nm, md->t, md->dt);
                }
                fds_exchange_divergence_info();
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

    if (useParallelPressure) {
        auto predPressureSubgraph = buildPressureIterationSubgraph(
            nmeshes, budget, true,
            depGraph, commService, fds_get_pres_flag());
        subgraph->input<TerminationData>(predPressureSubgraph);
        subgraph->edges(predDivP2KernelTask, predPressureSubgraph);
        subgraph->edges(predPressureSubgraph, velPredKernelTask);
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
    subgraph->edges(velPredKernelTask, changeTimeStepCollectorSM);
    subgraph->edges(changeTimeStepCollectorSM, changeTimeStepSubgraph);

    // Merged MeshExch(3) + SyntheticTurbulence -> PredFinal -> output
    subgraph->edges(changeTimeStepSubgraph, meshExch3SynTurbSM);
    subgraph->edges(meshExch3SynTurbSM, predFinalSubgraph);

    subgraph->outputs(predFinalSubgraph);

    return subgraph;
}

#endif // PREDICTOR_SUBGRAPH_H
