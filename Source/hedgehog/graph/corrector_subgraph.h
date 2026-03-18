#ifndef CORRECTOR_SUBGRAPH_H
#define CORRECTOR_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../state/collector_state.h"
#include "../state/div_setup_state.h"
#include "../state/velocity_corrector_state.h"
#include "../task/barrier_tasks.h"
#include "../task/corr_step1_kernel_task.h"
#include "../task/mass_fd_kernel_task.h"
#include "../task/div_setup_kernel_task.h"
#include "../task/combustion_kernel_task.h"
#include "../task/corr_condens_kernel_task.h"
#include "../task/particle_mass_energy_kernel_task.h"
#include "../task/corr_particle_kernel_task.h"
#include "../task/corr_div_part1_kernel_task.h"
#include "compute_viscosity_block_subgraph.h"
#include "velocity_flux_block_subgraph.h"
#include "../task/divergence_part2_kernel_task.h"
#include "divergence_part2_block_subgraph.h"
#include "../task/velocity_corrector_kernel_task.h"
#include "velocity_corrector_block_subgraph.h"
#include "particle_momentum_block_subgraph.h"
#include "wallbc_subgraph.h"
#include "wallbc_block_subgraph.h"
#include "velocity_bc_subgraph.h"
#include "corr_radiation_subgraph.h"
#include "pressure_iteration_subgraph.h"

/// Build the Corrector sub-graph.
///
/// Implements the full corrector phase of the FDS time-stepping loop:
///   CorrStep1 -> MESH_EXCHANGE(4) -> CorrDivSetup -> Combustion+HVAC ->
///   CorrCondens -> CorrParticle -> MESH_EXCHANGE(7) -> WallBC ->
///   MESH_EXCHANGE(6) -> CorrRadiation -> MESH_EXCHANGE(2)+InitDiv ->
///   CorrDivPart1 -> DivergenceExchange -> CorrDivPart2 -> PressureIteration ->
///   VelocityCorrector -> MESH_EXCHANGE(6) -> CorrFinal
///
/// Optimizations vs original graph:
///   - Combustion parallelized as kernel task, Soot+HVAC remains sequential barrier
///   - CorrRadiation outputs BarrierData directly (eliminates Collector(2))
///   - MeshExchange(2) includes InitDivIntegrals (eliminates 1 collector + 1 task)
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

    // --- Sub-graphs with orchestrators (sequential pre-processing required) ---

    // CorrDivSetup: parallel kernel (+ sequential CC_VELOCITY_BC if CC_IBM)
    // Block decomposition: if no Coriolis/patch/CTRL/wind/periodic, use K-block parallel
    bool ccIBM = fds_is_cc_ibm() != 0;
    bool canBlockFlux = fds_velocity_flux_can_block_decompose(1) != 0;
    auto corrDivSetupKernelTask = std::make_shared<DivSetupKernelTask>(kernelThreads);

    // CorrParticle: parallel MASS_ENERGY -> sequential REMOVE+MOVE -> parallel MOMENTUM
    auto particleMassEnergyKernelTask = std::make_shared<ParticleMassEnergyKernelTask>(kernelThreads);
    auto particleRemoveMoveCollSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "ParticleRemoveMoveCollector");
    auto particleRemoveMoveTask = std::make_shared<RemoveMoveParticlesTask>();
    // ParticleMomentum: block-decomposed for non-CC_IBM, mesh-level fallback for CC_IBM
    auto partMomSubgraph = buildParticleMomentumBlockSubgraph(
        blockThreads, numBlocks);
    auto corrParticleKernelTask = std::make_shared<CorrParticleKernelTask>(kernelThreads);

    // VelocityCorrector: block-decomposed kernel (+ CC_PROJECT_VELOCITY orch/collector if CC_IBM)
    // Block decomposition splits each mesh into K-range blocks for intra-mesh parallelism.
    // CHECK_DIVERGENCE_KERNEL runs at mesh level after block reassembly.
    // For sparse solvers (ULMAT/GLMAT/UGLMAT), WALL_VELOCITY_NO_GRADH runs inside
    // the block subgraph (store before decompose, fix after reassembly).
    auto velCorrSubgraph = buildVelocityCorrectorBlockSubgraph(
        blockThreads, numBlocks);
    // Fallback: original mesh-level kernel task for CC_IBM path
    auto velCorrKernelTask = std::make_shared<VelocityCorrectorKernelTask>(kernelThreads);

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

    // Combustion: parallel kernel -> Soot+HVAC sequential barrier
    auto combustionKernelTask = std::make_shared<CombustionKernelTask>(kernelThreads);
    auto sootHvacCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "SootHvacCollector");
    auto sootHvacTask = std::make_shared<SootHvacTask>(1);

    auto collector7SM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(7)");
    auto meshExchange7 = std::make_shared<MeshExchangeTask>(7);

    auto collector6aSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(6a)");
    auto meshExchange6a = std::make_shared<MeshExchangeTask>(6);

    // Merged: MeshExchange(2) + InitDivIntegrals (eliminates CorrInitDivCollector + InitDivTask)
    // CorrRadiation outputs BarrierData directly (eliminates Collector(2))
    auto meshExchange2 = std::make_shared<MeshExchangeTask>(2, /*ccDensity=*/false,
                                                             /*ccEndStep=*/false, /*initDiv=*/true);

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

    // CorrStep1: viscosity -> mass_fd + density -> MESH_EXCHANGE(4)
    if (canBlockVisc) {
        // Block-decomposed viscosity -> separate mass_fd + density task
        auto corrViscBlockSubgraph = buildComputeViscosityBlockSubgraph(
            nmeshes, blockThreads, numBlocks);
        auto corrMassFDDensityTask = std::make_shared<MassFDDensityKernelTask>(kernelThreads);
        subgraph->inputs(corrViscBlockSubgraph);
        subgraph->edges(corrViscBlockSubgraph, corrMassFDDensityTask);
        subgraph->edges(corrMassFDDensityTask, collector4SM);
    } else {
        // Mesh-level fallback: combined viscosity + mass_fd + density
        subgraph->inputs(corrStep1KernelTask);
        subgraph->edges(corrStep1KernelTask, collector4SM);
    }
    subgraph->edges(collector4SM, meshExchange4);

    // CorrDivSetup: block-decomposed or mesh-level depending on feature flags
    if (canBlockFlux) {
        // Block decomposition: orchestrator(pre-proc + K-decompose) -> parallel blocks -> collector
        // CC_IBM handled internally: CC_VELOCITY_BC + CUTFACE_VELOCITIES in orchestrator,
        // CC_VELOCITY_FLUX in collector
        auto corrDivSetupBlockSubgraph = buildVelocityFluxBlockSubgraph(
            nmeshes, blockThreads, numBlocks);
        subgraph->edges(meshExchange4, corrDivSetupBlockSubgraph);
        subgraph->edges(corrDivSetupBlockSubgraph, combustionKernelTask);
    } else if (ccIBM) {
        // CC_IBM with features preventing block decomposition: orchestrator + mesh-level kernel
        auto corrDivSetupOrchSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
            std::make_shared<CorrDivSetupOrchestrator>(nmeshes), "CorrDivSetupOrch");
        subgraph->edges(meshExchange4, corrDivSetupOrchSM);
        subgraph->edges(corrDivSetupOrchSM, corrDivSetupKernelTask);
        subgraph->edges(corrDivSetupKernelTask, combustionKernelTask);
    } else {
        // Mesh-level fallback (Coriolis, patch velocity, etc.)
        subgraph->edges(meshExchange4, corrDivSetupKernelTask);
        subgraph->edges(corrDivSetupKernelTask, combustionKernelTask);
    }

    // Combustion: parallel kernel -> Soot+HVAC barrier
    subgraph->edges(combustionKernelTask, sootHvacCollectorSM);
    subgraph->edges(sootHvacCollectorSM, sootHvacTask);

    // CorrCondens -> CorrParticle: parallel MASS_ENERGY -> REMOVE+MOVE barrier -> parallel MOMENTUM
    subgraph->edges(sootHvacTask, corrCondensKernelTask);
    subgraph->edges(corrCondensKernelTask, particleMassEnergyKernelTask);
    subgraph->edges(particleMassEnergyKernelTask, particleRemoveMoveCollSM);
    subgraph->edges(particleRemoveMoveCollSM, particleRemoveMoveTask);
    if (ccIBM) {
        subgraph->edges(particleRemoveMoveTask, corrParticleKernelTask);
        subgraph->edges(corrParticleKernelTask, collector7SM);
    } else {
        subgraph->edges(particleRemoveMoveTask, partMomSubgraph);
        subgraph->edges(partMomSubgraph, collector7SM);
    }
    subgraph->edges(collector7SM, meshExchange7);

    // WallBC sub-graph
    subgraph->edges(meshExchange7, wallBCSubgraph);
    subgraph->edges(wallBCSubgraph, collector6aSM);
    subgraph->edges(collector6aSM, meshExchange6a);

    // CorrRadiation sub-graph (outputs BarrierData directly)
    subgraph->edges(meshExchange6a, corrRadiationSubgraph);

    // MeshExchange(2) + InitDivIntegrals merged (CorrRadiation -> BarrierData -> MeshExchange2+InitDiv)
    subgraph->edges(corrRadiationSubgraph, meshExchange2);

    // CorrDivPart1 -> DivExchange
    subgraph->edges(meshExchange2, corrDivP1KernelTask);
    subgraph->edges(corrDivP1KernelTask, corrDivCollectorSM);
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
    // VelocityCorrector: parallel kernel (+ CC_PROJECT_VELOCITY orch/collector if CC_IBM)
    if (useParallelPressure) {
        auto corrPressureSubgraph = buildPressureIterationSubgraph(
            tEnd, nmeshes, kernelThreads, /*predictor=*/false, termSignal,
            fds_get_pres_flag());
        subgraph->edges(corrPressureCollectorSM, corrPressureSubgraph);
        // CC_IBM is always false when useParallelPressure is true
        subgraph->edges(corrPressureSubgraph, velCorrSubgraph);
        subgraph->edges(velCorrSubgraph, collector6bSM);
    } else {
        auto corrPressureTask = std::make_shared<PressureIterationTask>(/*predictor=*/false);
        subgraph->edges(corrPressureCollectorSM, corrPressureTask);
        if (ccIBM) {
            // CC_IBM uses mesh-level kernel (CC_PROJECT_VELOCITY requires full mesh)
            auto velCorrCCOrchSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
                std::make_shared<VelocityCorrectorCCOrchestrator>(nmeshes), "VelCorrCCOrch");
            auto velCorrCCCollSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
                std::make_shared<VelocityCorrectorCCCollector>(nmeshes), "VelCorrCCCollector");
            subgraph->edges(corrPressureTask, velCorrCCOrchSM);
            subgraph->edges(velCorrCCOrchSM, velCorrKernelTask);
            subgraph->edges(velCorrKernelTask, velCorrCCCollSM);
            subgraph->edges(velCorrCCCollSM, collector6bSM);
        } else {
            // Non-CC_IBM: use block-decomposed sub-graph for intra-mesh parallelism
            // (WallVelStore/Fix for sparse solvers handled inside the block subgraph)
            subgraph->edges(corrPressureTask, velCorrSubgraph);
            subgraph->edges(velCorrSubgraph, collector6bSM);
        }
    }
    subgraph->edges(collector6bSM, meshExchange6b);

    // CorrFinal sub-graph (outputs BarrierData directly — no external collector needed)
    subgraph->edges(meshExchange6b, corrFinalSubgraph);

    subgraph->outputs(corrFinalSubgraph);

    return subgraph;
}

#endif // CORRECTOR_SUBGRAPH_H
