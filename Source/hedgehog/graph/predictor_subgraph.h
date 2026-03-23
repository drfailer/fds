#ifndef PREDICTOR_SUBGRAPH_H
#define PREDICTOR_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../state/collector_state.h"
#include "../state/barrier_state.h"
#include "../state/pred_step1_state.h"
#include "../state/div_setup_state.h"
#include "../state/fork_join_state.h"
#include "../task/barrier_tasks.h"
#include "../task/pred_step1_kernel_task.h"
#include "../task/mass_fd_kernel_task.h"
#include "../task/density_pred_kernel_task.h"
#include "../task/div_setup_kernel_task.h"
#include "../task/pred_fork_tasks.h"
#include "../task/pred_wall_div_kernel_task.h"
#include "compute_viscosity_block_subgraph.h"
#include "velocity_flux_block_subgraph.h"
#include "../task/divergence_part2_kernel_task.h"
#include "divergence_part2_block_subgraph.h"
#include "velocity_predictor_block_subgraph.h"
#include "change_timestep_subgraph.h"
#include "velocity_bc_subgraph.h"
#include "wallbc_subgraph.h"
#include "wallbc_block_subgraph.h"
#include "pressure_iteration_subgraph.h"
#include "density_block_subgraph.h"
#include "visc_density_block_subgraph.h"
#include "pred_fork_vflux_subgraph.h"
#include "pred_fork_div_subgraph.h"

/// Build the Predictor sub-graph.
inline auto buildPredictorSubgraph(int nmeshes, double tEnd, size_t kernelThreads,
                                    size_t blockThreads, int numBlocks,
                                    std::shared_ptr<TerminationSignal> termSignal) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData, MeshData>>("Predictor");

    size_t meshThreads = static_cast<size_t>(nmeshes);

    // --- Kernel sub-graph components ---

    auto predStep1OrchSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<PredStep1Orchestrator>(nmeshes), "PredStep1Orch");
    auto predStep1KernelTask = std::make_shared<PredStep1KernelTask>(meshThreads);

    bool canBlockVisc = fds_compute_viscosity_can_block_decompose() != 0;
    bool canBlockDensity = fds_density_can_block_decompose() != 0 && numBlocks > 1;

    auto densPredKernelTask = std::make_shared<DensityPredKernelTask>(meshThreads);

    bool ccIBM = fds_is_cc_ibm() != 0;
    bool canBlockFlux = fds_velocity_flux_can_block_decompose(1) != 0;
    bool canBlockWallBC = fds_wall_bc_can_block_decompose() != 0;

    bool canBlockDivP2 = fds_divergence_part_2_can_block_decompose() != 0 && numBlocks > 1;
    auto predDivP2KernelTask = std::make_shared<DivergencePart2KernelTask>(meshThreads);

    auto velPredSubgraph = buildVelocityPredictorBlockSubgraph(
        blockThreads, numBlocks, nmeshes);

    auto predFinalSubgraph = buildPredFinalSubgraph(nmeshes, meshThreads, blockThreads, numBlocks);

    auto changeTimeStepSubgraph = buildChangeTimeStepSubgraph(tEnd, nmeshes, meshThreads);

    // --- Merged barrier states ---

    auto meshExchange1SM = makeBarrierSM(nmeshes, "MeshExchange(1)",
        "CC_DENSITY\\nMESH_EXCHANGE(1)\\nEXCHANGE_INSERTED_PARTICLES",
        [ccIBM](auto& meshes) {
            if (ccIBM) { fds_cc_density(meshes[0]->t, meshes[0]->dt); }
            fds_mesh_exchange(1);
            fds_exchange_inserted_particles();
        });

    auto predDivExchangeSM = makeBarrierSM(nmeshes, "PredDivExchange",
        "EXCHANGE_DIVERGENCE_INFO\\nGLOBAL_MATRIX_REASSIGN",
        [](auto& meshes) {
            fds_exchange_divergence_info();
            fds_global_matrix_reassign(0);
        });

    auto changeTimeStepCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "ChangeTimeStepCollector");

    auto meshExchange3SM = makeBarrierSM(nmeshes, "MeshExchange(3)",
        "CC_END_STEP\\nMESH_EXCHANGE(3)",
        [ccIBM](auto& meshes) {
            if (ccIBM) { fds_cc_end_step(meshes[0]->t, meshes[0]->dt, 0); }
            fds_mesh_exchange(3);
        });

    auto phaseTransTask = std::make_shared<PhaseTransitionTask>();

    // --- Wire the sub-graph ---

    subgraph->inputs(predStep1OrchSM);

    // PredStep1: orchestrator -> viscosity -> mass_fd -> Density
    if (canBlockVisc && canBlockDensity) {
        auto predViscDensitySubgraph = buildViscDensityBlockSubgraph(
            nmeshes, blockThreads, numBlocks);
        subgraph->edges(predStep1OrchSM, predViscDensitySubgraph);
        subgraph->edges(predViscDensitySubgraph, meshExchange1SM);
    } else if (canBlockVisc) {
        auto predViscBlockSubgraph = buildComputeViscosityBlockSubgraph(
            nmeshes, blockThreads, numBlocks);
        auto predMassFDKernelTask = std::make_shared<MassFDKernelTask>(meshThreads);
        subgraph->edges(predStep1OrchSM, predViscBlockSubgraph);
        subgraph->edges(predViscBlockSubgraph, predMassFDKernelTask);
        subgraph->edges(predMassFDKernelTask, densPredKernelTask);
        subgraph->edges(densPredKernelTask, meshExchange1SM);
    } else if (canBlockDensity) {
        auto predDensityBlockSubgraph = buildDensityBlockSubgraph(
            nmeshes, blockThreads, numBlocks);
        subgraph->edges(predStep1OrchSM, predStep1KernelTask);
        subgraph->edges(predStep1KernelTask, predDensityBlockSubgraph);
        subgraph->edges(predDensityBlockSubgraph, meshExchange1SM);
    } else {
        subgraph->edges(predStep1OrchSM, predStep1KernelTask);
        subgraph->edges(predStep1KernelTask, densPredKernelTask);
        subgraph->edges(densPredKernelTask, meshExchange1SM);
    }

    // --- Predictor middle section: Fork (non-CC_IBM) or Sequential (CC_IBM) ---

    if (!ccIBM) {
        // Merged barrier: HVAC + InitDiv + DivP1Prefork, then multicast to both branches
        auto hvacInitDivPreforkSM = makeBarrierSM(nmeshes, "HvacInitDivPrefork",
            "HVAC_CALC\\nINIT_DIV\\nDIV_P1_PREFORK",
            [](auto& meshes) {
                fds_hvac_calc(meshes[0]->t, meshes[0]->dt, 1);
                fds_initialize_divergence_integrals();
                for (auto &md : meshes) {
                    fds_divergence_part_1_prefork(md->nm, md->t, md->dt);
                }
            });

        subgraph->edges(meshExchange1SM, hvacInitDivPreforkSM);

        // Fork: (VFLUX + PART_MOM) || (WallBC + DIV_P1_early) via multicast
        auto predForkVFluxSG = buildPredForkVFluxSubgraph(
            nmeshes, blockThreads, numBlocks, canBlockFlux);
        auto predForkDivSG = buildPredForkDivSubgraph(
            nmeshes, blockThreads, numBlocks, canBlockWallBC);
        auto predJoinSM = std::make_shared<hh::StateManager<
            1, MeshData, MeshData>>(
            std::make_shared<ForkJoinState>(2), "PredJoin");

        // Multicast from barrier to both branches
        subgraph->edges(hvacInitDivPreforkSM, predForkVFluxSG);
        subgraph->edges(hvacInitDivPreforkSM, predForkDivSG);
        subgraph->edges(predForkVFluxSG, predJoinSM);
        subgraph->edges(predForkDivSG, predJoinSM);

        // After join: DIV_P1_late
        auto divP1LateTask = std::make_shared<DivP1LateTask>(meshThreads);
        subgraph->edges(predJoinSM, divP1LateTask);
        subgraph->edges(divP1LateTask, predDivExchangeSM);
    } else {
        // CC_IBM sequential path
        auto hvacInitDivSM = makeBarrierSM(nmeshes, "Hvac+InitDiv",
            "HVAC_CALC\\nINITIALIZE_DIVERGENCE_INTEGRALS",
            [](auto& meshes) {
                fds_hvac_calc(meshes[0]->t, meshes[0]->dt, 1);
                fds_initialize_divergence_integrals();
            });

        if (canBlockFlux) {
            auto predDivSetupBlockSubgraph = buildVelocityFluxBlockSubgraph(
                nmeshes, blockThreads, numBlocks);
            subgraph->edges(meshExchange1SM, predDivSetupBlockSubgraph);
            subgraph->edges(predDivSetupBlockSubgraph, hvacInitDivSM);
        } else {
            auto predDivSetupOrchSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
                std::make_shared<PredDivSetupOrchestrator>(nmeshes), "PredDivSetupOrch");
            auto predDivSetupKernelTask = std::make_shared<DivSetupKernelTask>(meshThreads);
            subgraph->edges(meshExchange1SM, predDivSetupOrchSM);
            subgraph->edges(predDivSetupOrchSM, predDivSetupKernelTask);
            subgraph->edges(predDivSetupKernelTask, hvacInitDivSM);
        }

        auto predWallBCSubgraph = canBlockWallBC
            ? buildWallBCBlockSubgraph(nmeshes, blockThreads, numBlocks)
            : buildWallBCSubgraph(nmeshes, meshThreads);
        subgraph->edges(hvacInitDivSM, predWallBCSubgraph);

        auto predWallDivKernelTask = std::make_shared<PredWallDivKernelTask>(meshThreads);
        subgraph->edges(predWallBCSubgraph, predWallDivKernelTask);
        subgraph->edges(predWallDivKernelTask, predDivExchangeSM);
    }

    // --- Common downstream: DivExchange -> DivP2 -> Pressure -> VelPred -> ... ---

    bool useParallelPressure = fds_use_pressure_subgraph() != 0;
    auto wirePressure = [&](auto lastDivP2Node) {
        if (useParallelPressure) {
            auto predPressureCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
                std::make_shared<CollectorState>(nmeshes), "PredPressureCollector");
            auto predPressureSubgraph = buildPressureIterationSubgraph(
                tEnd, nmeshes, meshThreads, true, termSignal,
                fds_get_pres_flag());
            subgraph->edges(lastDivP2Node, predPressureCollectorSM);
            subgraph->edges(predPressureCollectorSM, predPressureSubgraph);
            subgraph->edges(predPressureSubgraph, velPredSubgraph);
        } else {
            auto predPressureSM = makeBarrierSM(nmeshes, "PredPressure",
                "PRESSURE_ITERATION\\nINIT_CHANGE_TIME_STEP",
                [](auto& meshes) {
                    fds_pressure_iteration(meshes[0]->t, meshes[0]->dt);
                    fds_init_change_time_step(meshes[0]->dt);
                });
            subgraph->edges(lastDivP2Node, predPressureSM);
            subgraph->edges(predPressureSM, velPredSubgraph);
        }
    };

    if (canBlockDivP2) {
        auto predDivP2BlockSubgraph = buildDivergencePart2BlockSubgraph(
            nmeshes, blockThreads, numBlocks);
        subgraph->edges(predDivExchangeSM, predDivP2BlockSubgraph);
        wirePressure(predDivP2BlockSubgraph);
    } else {
        subgraph->edges(predDivExchangeSM, predDivP2KernelTask);
        wirePressure(predDivP2KernelTask);
    }

    // VelocityPredictor -> ChangeTimeStep
    subgraph->edges(velPredSubgraph, changeTimeStepCollectorSM);
    subgraph->edges(changeTimeStepCollectorSM, changeTimeStepSubgraph);
    subgraph->edges(changeTimeStepSubgraph, meshExchange3SM);

    // PredFinal -> PhaseTransition
    subgraph->edges(meshExchange3SM, predFinalSubgraph);
    subgraph->edges(predFinalSubgraph, phaseTransTask);

    subgraph->outputs(phaseTransTask);

    return subgraph;
}

#endif // PREDICTOR_SUBGRAPH_H
