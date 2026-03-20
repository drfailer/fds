#ifndef CORRECTOR_SUBGRAPH_H
#define CORRECTOR_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../state/collector_state.h"
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
#include "../data/pipeline_fork2_data.h"
#include "../state/pipeline_fork2_state.h"
#include "../task/pipeline_fork2_tasks.h"
#include "pipeline_fork2_rad_subgraph.h"

/// Build the Corrector sub-graph.
///
/// Implements the full corrector phase of the FDS time-stepping loop:
///   CorrStep1 -> MESH_EXCHANGE(4) -> CorrDivSetup -> Combustion+HVAC ->
///   CorrCondens -> CorrParticle -> MESH_EXCHANGE(7) -> WallBC ->
///   MESH_EXCHANGE(6) -> Fork2(RADIATION || DIV_P1_noQR) -> MESH_EXCHANGE(2) ->
///   QR_Addition -> DivergenceExchange -> CorrDivPart2 -> PressureIteration ->
///   VelocityCorrector -> MESH_EXCHANGE(6) -> CorrFinal
///
/// Optimizations vs original graph:
///   - Combustion parallelized as kernel task, Soot+HVAC remains sequential barrier
///   - CorrRadiation outputs BarrierData directly (eliminates Collector(2))
///   - Fork 2: RADIATION || DIV_P1(SKIP_QR, WORK_BRANCH=2) concurrent execution
///   - InitDivIntegrals in Fork2State; QR addition after MeshExchange(2)
///   - CC_IBM falls back to sequential RADIATION -> DIV_P1 (WORK_BRANCH=2 incompatible)
///   - CorrFinal outputs BarrierData directly (eliminates TimestepCollector in parent)
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

    // CorrDivSetup: parallel kernel (+ sequential CC_VELOCITY_BC if CC_IBM)
    // Block decomposition: if no Coriolis/patch/CTRL/wind/periodic, use K-block parallel
    bool ccIBM = fds_is_cc_ibm() != 0;
    bool canBlockFlux = fds_velocity_flux_can_block_decompose(1) != 0;
    // CorrParticle: parallel MASS_ENERGY -> sequential REMOVE+MOVE -> parallel MOMENTUM
    auto particleMassEnergyKernelTask = std::make_shared<ParticleMassEnergyKernelTask>(kernelThreads);
    auto particleRemoveMoveCollSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "ParticleRemoveMoveCollector");
    auto particleRemoveMoveTask = std::make_shared<RemoveMoveParticlesTask>();
    // ParticleMomentum: always block-decomposed (K-blocks are particle-safe)
    auto partMomSubgraph = buildParticleMomentumBlockSubgraph(
        blockThreads, numBlocks);

    // VelocityCorrector: block-decomposed kernel
    // CC_PROJECT_VELOCITY and WALL_VELOCITY_NO_GRADH are no-op tasks in the pipeline
    // for non-CC_IBM / FFT respectively (checked in Fortran C wrapper).
    // CHECK_DIVERGENCE_KERNEL runs at mesh level after block reassembly.
    auto velCorrSubgraph = buildVelocityCorrectorBlockSubgraph(
        blockThreads, numBlocks);

    // --- Named sub-graphs ---

    bool canBlockWallBC = fds_wall_bc_can_block_decompose() != 0;
    auto wallBCSubgraph = canBlockWallBC
        ? buildWallBCBlockSubgraph(nmeshes, blockThreads, numBlocks)
        : buildWallBCSubgraph(nmeshes, kernelThreads);
    auto corrRadiationSubgraph = buildCorrRadiationSubgraph(nmeshes, kernelThreads);
    auto corrFinalSubgraph = buildCorrFinalSubgraph(nmeshes, kernelThreads, blockThreads, numBlocks);

    // --- Barrier tasks ---

    auto collector4SM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(4)");
    auto meshExchange4 = std::make_shared<MeshExchangeTask>(4, /*ccDensity=*/ccIBM);

    auto sootHvacCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "SootHvacCollector");
    auto sootHvacTask = std::make_shared<SootHvacTask>(1);

    auto collector7SM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(7)");
    auto meshExchange7 = std::make_shared<MeshExchangeTask>(7);

    auto collector6aSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(6a)");
    auto meshExchange6a = std::make_shared<MeshExchangeTask>(6);

    // MeshExchange(2): QR exchange after radiation.
    // CC_IBM: sequential path keeps InitDiv in MeshExchange(2).
    // Non-CC_IBM: Fork2State handles InitDiv before dispatching branches.
    auto meshExchange2 = std::make_shared<MeshExchangeTask>(2, /*ccDensity=*/false,
                                                             /*ccEndStep=*/false, /*initDiv=*/ccIBM);

    auto corrDivCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "CorrDivCollector");
    auto corrDivExchangeTask = std::make_shared<DivergenceExchangeTask>(/*corrector=*/true);

    auto corrPressureCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "CorrPressureCollector");
    bool useParallelPressure = fds_use_pressure_subgraph() != 0;

    auto collector6bSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(6b)");
    auto meshExchange6b = std::make_shared<MeshExchangeTask>(6, /*ccDensity=*/false, /*ccEndStep=*/ccIBM);

    // --- Wire the sub-graph ---

    // CorrStep1: viscosity -> mass_fd -> density -> MESH_EXCHANGE(4)
    if (canBlockVisc) {
        // Block-decomposed viscosity -> separate mass_fd task
        auto corrViscBlockSubgraph = buildComputeViscosityBlockSubgraph(
            nmeshes, blockThreads, numBlocks);
        auto corrMassFDKernelTask = std::make_shared<MassFDKernelTask>(kernelThreads);
        subgraph->inputs(corrViscBlockSubgraph);
        subgraph->edges(corrViscBlockSubgraph, corrMassFDKernelTask);
        if (canBlockDensity) {
            auto corrDensityBlockSubgraph = buildDensityBlockSubgraph(
                nmeshes, blockThreads, numBlocks);
            subgraph->edges(corrMassFDKernelTask, corrDensityBlockSubgraph);
            subgraph->edges(corrDensityBlockSubgraph, collector4SM);
        } else {
            auto corrMassFDDensityFallback = std::make_shared<DensityPredKernelTask>(kernelThreads);
            subgraph->edges(corrMassFDKernelTask, corrMassFDDensityFallback);
            subgraph->edges(corrMassFDDensityFallback, collector4SM);
        }
    } else {
        if (canBlockDensity) {
            // Mesh-level viscosity + mass_fd, block-decomposed density
            auto corrViscMassFDTask = std::make_shared<CorrViscMassFDKernelTask>(kernelThreads);
            auto corrDensityBlockSubgraph = buildDensityBlockSubgraph(
                nmeshes, blockThreads, numBlocks);
            subgraph->inputs(corrViscMassFDTask);
            subgraph->edges(corrViscMassFDTask, corrDensityBlockSubgraph);
            subgraph->edges(corrDensityBlockSubgraph, collector4SM);
        } else {
            // Mesh-level fallback: combined viscosity + mass_fd + density
            subgraph->inputs(corrStep1KernelTask);
            subgraph->edges(corrStep1KernelTask, collector4SM);
        }
    }
    subgraph->edges(collector4SM, meshExchange4);

    // --- Fork 1: VFLUX || COMBUSTION ---
    // Both branches run concurrently after MeshExchange(4), join before SootHvac.

    auto fork1SM = std::make_shared<hh::StateManager<1, MeshData, Fork1VFluxWork, Fork1CombWork>>(
        std::make_shared<PipelineFork1State>(), "Fork1");
    auto fork1VFluxSubgraph = buildFork1VFluxSubgraph(
        nmeshes, kernelThreads, blockThreads, numBlocks, canBlockFlux, ccIBM);
    auto fork1CombTask = std::make_shared<Fork1CombKernelTask>(kernelThreads);
    auto join1SM = std::make_shared<hh::StateManager<2, Fork1VFluxResult, Fork1CombResult, MeshData>>(
        std::make_shared<PipelineJoin1State>(), "Join1");

    subgraph->edges(meshExchange4, fork1SM);
    // Branch A: Fork1 -> VFLUX sub-graph -> Join1
    subgraph->edges(fork1SM, fork1VFluxSubgraph);
    subgraph->edges(fork1VFluxSubgraph, join1SM);
    // Branch B: Fork1 -> Combustion task -> Join1
    subgraph->edges(fork1SM, fork1CombTask);
    subgraph->edges(fork1CombTask, join1SM);

    // After join: Soot+HVAC barrier
    subgraph->edges(join1SM, sootHvacCollectorSM);
    subgraph->edges(sootHvacCollectorSM, sootHvacTask);

    // CorrCondens -> CorrParticle: parallel MASS_ENERGY -> REMOVE+MOVE barrier -> parallel MOMENTUM
    subgraph->edges(sootHvacTask, corrCondensKernelTask);
    subgraph->edges(corrCondensKernelTask, particleMassEnergyKernelTask);
    subgraph->edges(particleMassEnergyKernelTask, particleRemoveMoveCollSM);
    subgraph->edges(particleRemoveMoveCollSM, particleRemoveMoveTask);
    subgraph->edges(particleRemoveMoveTask, partMomSubgraph);
    subgraph->edges(partMomSubgraph, collector7SM);
    subgraph->edges(collector7SM, meshExchange7);

    // WallBC sub-graph
    subgraph->edges(meshExchange7, wallBCSubgraph);
    subgraph->edges(wallBCSubgraph, collector6aSM);
    subgraph->edges(collector6aSM, meshExchange6a);

    // --- Fork 2: RADIATION || DIV_P1 (or sequential fallback for CC_IBM) ---
    if (ccIBM) {
        // CC_IBM: sequential path (WORK_BRANCH=2 incompatible with CC divergence code)
        subgraph->edges(meshExchange6a, corrRadiationSubgraph);
        subgraph->edges(corrRadiationSubgraph, meshExchange2);
        subgraph->edges(meshExchange2, corrDivP1KernelTask);
        subgraph->edges(corrDivP1KernelTask, corrDivCollectorSM);
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

        subgraph->edges(meshExchange6a, fork2SM);
        // Branch C: Radiation
        subgraph->edges(fork2SM, fork2RadSubgraph);
        subgraph->edges(fork2RadSubgraph, join2SM);
        // Branch D: DIV_P1 (SKIP_QR, WORK_BRANCH=2)
        subgraph->edges(fork2SM, fork2DivP1Task);
        subgraph->edges(fork2DivP1Task, fork2DivP1CollSM);
        subgraph->edges(fork2DivP1CollSM, join2SM);
        // After join: MeshExchange(2) -> QR addition -> collector
        subgraph->edges(join2SM, meshExchange2);
        subgraph->edges(meshExchange2, qrAddTask);
        subgraph->edges(qrAddTask, corrDivCollectorSM);
    }
    subgraph->edges(corrDivCollectorSM, corrDivExchangeTask);

    // CorrDivPart2 -> Pressure (block-decomposed or mesh-level)
    if (canBlockDivP2) {
        auto corrDivP2BlockSubgraph = buildDivergencePart2BlockSubgraph(
            nmeshes, blockThreads, numBlocks);
        subgraph->edges(corrDivExchangeTask, corrDivP2BlockSubgraph);
        subgraph->edges(corrDivP2BlockSubgraph, corrPressureCollectorSM);
    } else {
        subgraph->edges(corrDivExchangeTask, corrDivP2KernelTask);
        subgraph->edges(corrDivP2KernelTask, corrPressureCollectorSM);
    }

    // Pressure iteration: parallel sub-graph or sequential fallback
    if (useParallelPressure) {
        auto corrPressureSubgraph = buildPressureIterationSubgraph(
            tEnd, nmeshes, kernelThreads, /*predictor=*/false, termSignal,
            fds_get_pres_flag());
        subgraph->edges(corrPressureCollectorSM, corrPressureSubgraph);
        subgraph->edges(corrPressureSubgraph, velCorrSubgraph);
    } else {
        auto corrPressureTask = std::make_shared<PressureIterationTask>(/*predictor=*/false);
        subgraph->edges(corrPressureCollectorSM, corrPressureTask);
        subgraph->edges(corrPressureTask, velCorrSubgraph);
    }
    subgraph->edges(velCorrSubgraph, collector6bSM);
    subgraph->edges(collector6bSM, meshExchange6b);

    // CorrFinal sub-graph (outputs BarrierData directly — no external collector needed)
    subgraph->edges(meshExchange6b, corrFinalSubgraph);

    subgraph->outputs(corrFinalSubgraph);

    return subgraph;
}

#endif // CORRECTOR_SUBGRAPH_H
