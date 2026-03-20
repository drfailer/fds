#ifndef CORRECTOR_SUBGRAPH_H
#define CORRECTOR_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../state/collector_state.h"
#include "../state/barrier_state.h"
#include "../state/div_setup_state.h"
#include "../task/barrier_tasks.h"
#include "../task/corr_step1_kernel_task.h"
#include "../task/mass_fd_kernel_task.h"
#include "../task/div_setup_kernel_task.h"
#include "../task/combustion_kernel_task.h"
#include "../data/pipeline_fork1_data.h"
#include "../state/pipeline_fork1_state.h"
#include "../task/pipeline_fork1_tasks.h"
#include "pipeline_fork1_vflux_subgraph.h"
#include "../task/corr_condens_kernel_task.h"
#include "../task/particle_mass_energy_kernel_task.h"
#include "../task/corr_div_part1_kernel_task.h"
#include "compute_viscosity_block_subgraph.h"
#include "velocity_flux_block_subgraph.h"
#include "../task/divergence_part2_kernel_task.h"
#include "divergence_part2_block_subgraph.h"
#include "velocity_corrector_block_subgraph.h"
#include "particle_momentum_block_subgraph.h"
#include "wallbc_subgraph.h"
#include "wallbc_block_subgraph.h"
#include "velocity_bc_subgraph.h"
#include "corr_radiation_subgraph.h"
#include "pressure_iteration_subgraph.h"
#include "density_block_subgraph.h"
#include "visc_density_block_subgraph.h"
#include "../data/pipeline_fork2_data.h"
#include "../state/pipeline_fork2_state.h"
#include "../task/pipeline_fork2_tasks.h"
#include "pipeline_fork2_rad_subgraph.h"

/// Build the Corrector sub-graph.
///
/// Implements the full corrector phase of the FDS time-stepping loop:
///   CorrStep1 -> MESH_EXCHANGE(4) -> Fork1(VFLUX || COMB) -> Soot+HVAC ->
///   CorrCondens -> CorrParticle -> MESH_EXCHANGE(7) -> WallBC ->
///   MESH_EXCHANGE(6a) -> Fork2(RADIATION || DIV_P1_noQR) -> MESH_EXCHANGE(2) ->
///   QR_Addition -> DivExchange -> CorrDivPart2 -> PressureIteration ->
///   VelocityCorrector -> MESH_EXCHANGE(6b) -> CorrFinal
///
/// @param nmeshes Number of meshes
/// @param tEnd Simulation end time (for pressure iteration sub-graph termination)
/// @param kernelThreads Number of threads for parallel kernel tasks
/// @param termSignal Shared termination signal for pressure iteration sub-graph
/// @return Shared pointer to the constructed sub-graph
inline auto buildCorrectorSubgraph(int nmeshes, double tEnd, size_t kernelThreads,
                                    size_t blockThreads, int numBlocks,
                                    std::shared_ptr<TerminationSignal> termSignal) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData, BarrierData>>("Corrector");

    // --- Kernel tasks (MeshData -> MeshData, no orchestrator/collector needed) ---

    auto corrStep1KernelTask = std::make_shared<CorrStep1KernelTask>(kernelThreads);
    auto corrCondensKernelTask = std::make_shared<CorrCondensKernelTask>(kernelThreads);
    auto corrDivP1KernelTask = std::make_shared<CorrDivPart1KernelTask>(kernelThreads);
    auto corrDivP2KernelTask = std::make_shared<DivergencePart2KernelTask>(kernelThreads);
    bool canBlockDivP2 = fds_divergence_part_2_can_block_decompose() != 0 && numBlocks > 1;

    // Viscosity block decomposition: if non-DEARDORFF/DYNSMAG/CC_IBM, use K-block parallel
    bool canBlockVisc = fds_compute_viscosity_can_block_decompose() != 0;

    // Density block decomposition: if non-CC_IBM and non-MMS, use K-block parallel
    bool canBlockDensity = fds_density_can_block_decompose() != 0 && numBlocks > 1;

    // --- Sub-graphs with orchestrators (sequential pre-processing required) ---

    bool ccIBM = fds_is_cc_ibm() != 0;
    bool canBlockFlux = fds_velocity_flux_can_block_decompose(1) != 0;

    // CorrParticle: parallel MASS_ENERGY -> sequential REMOVE+MOVE -> parallel MOMENTUM
    auto particleMassEnergyKernelTask = std::make_shared<ParticleMassEnergyKernelTask>(kernelThreads);
    auto partMomSubgraph = buildParticleMomentumBlockSubgraph(
        blockThreads, numBlocks);

    // VelocityCorrector: block-decomposed kernel
    auto velCorrSubgraph = buildVelocityCorrectorBlockSubgraph(
        blockThreads, numBlocks);

    // --- Named sub-graphs ---

    bool canBlockWallBC = fds_wall_bc_can_block_decompose() != 0;
    auto wallBCSubgraph = canBlockWallBC
        ? buildWallBCBlockSubgraph(nmeshes, blockThreads, numBlocks)
        : buildWallBCSubgraph(nmeshes, kernelThreads);
    auto corrRadiationSubgraph = buildCorrRadiationSubgraph(nmeshes, kernelThreads);
    auto corrFinalSubgraph = buildCorrFinalSubgraph(nmeshes, kernelThreads, blockThreads, numBlocks);

    // --- Merged barrier states (replace collector + barrier task pairs) ---

    auto meshExchange4SM = makeBarrierSM(nmeshes, "MeshExchange(4)",
        "CC_DENSITY\\nMESH_EXCHANGE(4)",
        [ccIBM](auto& meshes) {
            if (ccIBM) { fds_cc_density(meshes[0]->t, meshes[0]->dt); }
            fds_mesh_exchange(4);
        });

    auto sootHvacSM = makeBarrierSM(nmeshes, "Soot+Hvac",
        "SOOT_OXIDATION_LOOP\\nHVAC_CALC",
        [](auto& meshes) {
            fds_soot_oxidation_loop(meshes[0]->dt);
            fds_hvac_calc(meshes[0]->t, meshes[0]->dt, 1);
        });

    auto removeMoveSM = makeBarrierSM(nmeshes, "RemoveMove",
        "REMOVE_PARTICLES\\nMOVE_PARTICLES",
        [](auto& meshes) {
            for (auto &md : meshes) {
                fds_remove_particles(md->t, md->nm);
                fds_move_particles(md->t, md->dt, md->nm);
            }
        });

    auto meshExchange7SM = makeBarrierSM(nmeshes, "MeshExchange(7)",
        "MESH_EXCHANGE(7)",
        [](auto& meshes) { fds_mesh_exchange(7); });

    auto meshExchange6aSM = makeBarrierSM(nmeshes, "MeshExchange(6a)",
        "MESH_EXCHANGE(6)",
        [](auto& meshes) { fds_mesh_exchange(6); });

    // MeshExchange(2): standalone BarrierData->MeshData task (no collector to merge).
    // CC_IBM: sequential path keeps InitDiv in MeshExchange(2).
    // Non-CC_IBM: Fork2State handles InitDiv before dispatching branches.
    auto meshExchange2 = std::make_shared<MeshExchangeTask>(2, /*ccDensity=*/false,
                                                             /*ccEndStep=*/false, /*initDiv=*/ccIBM);

    auto corrDivExchangeSM = makeBarrierSM(nmeshes, "CorrDivExchange",
        "EXCHANGE_DIVERGENCE_INFO\\nRTE_SOURCE_CORRECTION\\nGLOBAL_MATRIX_REASSIGN",
        [](auto& meshes) {
            fds_exchange_divergence_info();
            fds_rte_source_correction();
            fds_global_matrix_reassign(0);
        });

    bool useParallelPressure = fds_use_pressure_subgraph() != 0;

    auto meshExchange6bSM = makeBarrierSM(nmeshes, "MeshExchange(6b)",
        "CC_END_STEP\\nMESH_EXCHANGE(6)",
        [ccIBM](auto& meshes) {
            if (ccIBM) { fds_cc_end_step(meshes[0]->t, meshes[0]->dt, 0); }
            fds_mesh_exchange(6);
        });

    // --- Wire the sub-graph ---

    // CorrStep1: viscosity -> mass_fd -> density -> MESH_EXCHANGE(4)
    if (canBlockVisc && canBlockDensity) {
        // Merged: ViscBlockKernel -> MidState(ViscPost+MassFD+DensPrep) -> DensityBlockKernel
        auto corrViscDensitySubgraph = buildViscDensityBlockSubgraph(
            nmeshes, blockThreads, numBlocks);
        subgraph->inputs(corrViscDensitySubgraph);
        subgraph->edges(corrViscDensitySubgraph, meshExchange4SM);
    } else if (canBlockVisc) {
        auto corrViscBlockSubgraph = buildComputeViscosityBlockSubgraph(
            nmeshes, blockThreads, numBlocks);
        auto corrMassFDKernelTask = std::make_shared<MassFDKernelTask>(kernelThreads);
        auto corrMassFDDensityFallback = std::make_shared<DensityPredKernelTask>(kernelThreads);
        subgraph->inputs(corrViscBlockSubgraph);
        subgraph->edges(corrViscBlockSubgraph, corrMassFDKernelTask);
        subgraph->edges(corrMassFDKernelTask, corrMassFDDensityFallback);
        subgraph->edges(corrMassFDDensityFallback, meshExchange4SM);
    } else if (canBlockDensity) {
        auto corrViscMassFDTask = std::make_shared<CorrViscMassFDKernelTask>(kernelThreads);
        auto corrDensityBlockSubgraph = buildDensityBlockSubgraph(
            nmeshes, blockThreads, numBlocks);
        subgraph->inputs(corrViscMassFDTask);
        subgraph->edges(corrViscMassFDTask, corrDensityBlockSubgraph);
        subgraph->edges(corrDensityBlockSubgraph, meshExchange4SM);
    } else {
        subgraph->inputs(corrStep1KernelTask);
        subgraph->edges(corrStep1KernelTask, meshExchange4SM);
    }

    // --- Fork 1: VFLUX || COMBUSTION ---

    auto fork1SM = std::make_shared<hh::StateManager<1, MeshData, Fork1VFluxWork, Fork1CombWork>>(
        std::make_shared<PipelineFork1State>(), "Fork1");
    auto fork1VFluxSubgraph = buildFork1VFluxSubgraph(
        nmeshes, kernelThreads, blockThreads, numBlocks, canBlockFlux, ccIBM);
    auto fork1CombTask = std::make_shared<Fork1CombKernelTask>(kernelThreads);
    auto join1SM = std::make_shared<hh::StateManager<2, Fork1VFluxResult, Fork1CombResult, MeshData>>(
        std::make_shared<PipelineJoin1State>(), "Join1");

    subgraph->edges(meshExchange4SM, fork1SM);
    // Branch A: Fork1 -> VFLUX sub-graph -> Join1
    subgraph->edges(fork1SM, fork1VFluxSubgraph);
    subgraph->edges(fork1VFluxSubgraph, join1SM);
    // Branch B: Fork1 -> Combustion task -> Join1
    subgraph->edges(fork1SM, fork1CombTask);
    subgraph->edges(fork1CombTask, join1SM);

    // After join: Soot+HVAC barrier
    subgraph->edges(join1SM, sootHvacSM);

    // CorrCondens -> CorrParticle: parallel MASS_ENERGY -> REMOVE+MOVE barrier -> parallel MOMENTUM
    subgraph->edges(sootHvacSM, corrCondensKernelTask);
    subgraph->edges(corrCondensKernelTask, particleMassEnergyKernelTask);
    subgraph->edges(particleMassEnergyKernelTask, removeMoveSM);
    subgraph->edges(removeMoveSM, partMomSubgraph);
    subgraph->edges(partMomSubgraph, meshExchange7SM);

    // WallBC sub-graph
    subgraph->edges(meshExchange7SM, wallBCSubgraph);
    subgraph->edges(wallBCSubgraph, meshExchange6aSM);

    // --- Fork 2: RADIATION || DIV_P1 (or sequential fallback for CC_IBM) ---
    if (ccIBM) {
        // CC_IBM: sequential path (WORK_BRANCH=2 incompatible with CC divergence code)
        subgraph->edges(meshExchange6aSM, corrRadiationSubgraph);
        subgraph->edges(corrRadiationSubgraph, meshExchange2);
        subgraph->edges(meshExchange2, corrDivP1KernelTask);
        subgraph->edges(corrDivP1KernelTask, corrDivExchangeSM);
    } else {
        // Fork 2: RADIATION (WORK_BRANCH=1) || DIV_P1 (SKIP_QR, WORK_BRANCH=2)
        auto fork2SM = std::make_shared<hh::StateManager<
            1, MeshData, Fork2RadWork, Fork2DivP1Work>>(
            std::make_shared<PipelineFork2State>(nmeshes), "Fork2");
        auto fork2RadSubgraph = buildFork2RadSubgraph(nmeshes, kernelThreads);
        auto fork2DivP1Task = std::make_shared<Fork2DivP1KernelTask>(kernelThreads);
        auto fork2DivP1CollSM = std::make_shared<hh::StateManager<
            1, Fork2DivP1Work, Fork2DivP1Barrier>>(
            std::make_shared<Fork2DivP1CollectorState>(nmeshes),
            "Fork2DivP1Collector");
        auto join2SM = std::make_shared<hh::StateManager<
            2, Fork2RadBarrier, Fork2DivP1Barrier, BarrierData>>(
            std::make_shared<PipelineJoin2State>(), "Join2");
        auto qrAddTask = std::make_shared<DivP1QRAdditionTask>(kernelThreads);

        subgraph->edges(meshExchange6aSM, fork2SM);
        // Branch C: Radiation
        subgraph->edges(fork2SM, fork2RadSubgraph);
        subgraph->edges(fork2RadSubgraph, join2SM);
        // Branch D: DIV_P1 (SKIP_QR, WORK_BRANCH=2)
        subgraph->edges(fork2SM, fork2DivP1Task);
        subgraph->edges(fork2DivP1Task, fork2DivP1CollSM);
        subgraph->edges(fork2DivP1CollSM, join2SM);
        // After join: MeshExchange(2) -> QR addition -> DivExchange
        subgraph->edges(join2SM, meshExchange2);
        subgraph->edges(meshExchange2, qrAddTask);
        subgraph->edges(qrAddTask, corrDivExchangeSM);
    }

    // --- Common downstream: DivExchange -> DivP2 -> Pressure -> VelCorr -> ... ---

    // Helper lambda to wire from DivP2 output through pressure to velCorrSubgraph
    auto wirePressure = [&](auto lastDivP2Node) {
        if (useParallelPressure) {
            auto corrPressureCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
                std::make_shared<CollectorState>(nmeshes), "CorrPressureCollector");
            auto corrPressureSubgraph = buildPressureIterationSubgraph(
                tEnd, nmeshes, kernelThreads, /*predictor=*/false, termSignal,
                fds_get_pres_flag());
            subgraph->edges(lastDivP2Node, corrPressureCollectorSM);
            subgraph->edges(corrPressureCollectorSM, corrPressureSubgraph);
            subgraph->edges(corrPressureSubgraph, velCorrSubgraph);
        } else {
            auto corrPressureSM = makeBarrierSM(nmeshes, "CorrPressure",
                "PRESSURE_ITERATION",
                [](auto& meshes) {
                    fds_pressure_iteration(meshes[0]->t, meshes[0]->dt);
                });
            subgraph->edges(lastDivP2Node, corrPressureSM);
            subgraph->edges(corrPressureSM, velCorrSubgraph);
        }
    };

    if (canBlockDivP2) {
        auto corrDivP2BlockSubgraph = buildDivergencePart2BlockSubgraph(
            nmeshes, blockThreads, numBlocks);
        subgraph->edges(corrDivExchangeSM, corrDivP2BlockSubgraph);
        wirePressure(corrDivP2BlockSubgraph);
    } else {
        subgraph->edges(corrDivExchangeSM, corrDivP2KernelTask);
        wirePressure(corrDivP2KernelTask);
    }

    subgraph->edges(velCorrSubgraph, meshExchange6bSM);

    // CorrFinal sub-graph (outputs BarrierData directly)
    subgraph->edges(meshExchange6bSM, corrFinalSubgraph);

    subgraph->outputs(corrFinalSubgraph);

    return subgraph;
}

#endif // CORRECTOR_SUBGRAPH_H
