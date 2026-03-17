#ifndef PREDICTOR_SUBGRAPH_H
#define PREDICTOR_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../state/collector_state.h"
#include "../state/pred_step1_state.h"
#include "../state/div_setup_state.h"
#include "../state/velocity_predictor_state.h"
#include "../task/barrier_tasks.h"
#include "../task/pred_step1_kernel_task.h"
#include "../task/mass_fd_kernel_task.h"
#include "../task/density_pred_kernel_task.h"
#include "../task/div_setup_kernel_task.h"
#include "../task/pred_wall_div_kernel_task.h"
#include "compute_viscosity_block_subgraph.h"
#include "velocity_flux_block_subgraph.h"
#include "../task/divergence_part2_kernel_task.h"
#include "divergence_part2_block_subgraph.h"
#include "../task/velocity_predictor_kernel_task.h"
#include "velocity_predictor_block_subgraph.h"
#include "change_timestep_subgraph.h"
#include "velocity_bc_subgraph.h"
#include "wallbc_subgraph.h"
#include "wallbc_block_subgraph.h"
#include "pressure_iteration_subgraph.h"

/// Build the Predictor sub-graph.
///
/// Implements the full predictor phase of the FDS time-stepping loop:
///   PredStep1 -> DensityPred -> MESH_EXCHANGE(1) -> PredDivSetup -> HVAC+InitDiv ->
///   WallBC -> PredWallDiv -> DivergenceExchange -> PredDivPart2 ->
///   PressureIteration -> VelocityPredictor -> ChangeTimeStep ->
///   MESH_EXCHANGE(3) -> PredFinal -> PhaseTransition
///
/// Optimizations vs original graph:
///   - HVAC + InitDivIntegrals merged into single barrier task (eliminates 1 collector + 1 task)
///   - PredFinal outputs BarrierData directly (eliminates PhaseTransCollector)
///
/// @param nmeshes Number of meshes
/// @param tEnd Simulation end time (passed to ChangeTimeStep sub-graph)
/// @param kernelThreads Number of threads for parallel kernel tasks
/// @param termSignal Shared termination signal for pressure iteration sub-graph
/// @return Shared pointer to the constructed sub-graph
inline auto buildPredictorSubgraph(int nmeshes, double tEnd, size_t kernelThreads,
                                    size_t blockThreads, int numBlocks,
                                    std::shared_ptr<TerminationSignal> termSignal) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData, MeshData>>("Predictor");

    // --- Kernel sub-graph components ---

    // PredStep1: sequential INSERT_ALL_PARTICLES -> parallel kernels
    auto predStep1OrchSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<PredStep1Orchestrator>(nmeshes), "PredStep1Orch");
    auto predStep1KernelTask = std::make_shared<PredStep1KernelTask>(kernelThreads);

    // Viscosity block decomposition: if non-DEARDORFF/DYNSMAG/CC_IBM, use K-block parallel
    bool canBlockVisc = fds_compute_viscosity_can_block_decompose() != 0;

    // DensityPred: parallel DENSITY_KERNEL
    auto densPredKernelTask = std::make_shared<DensityPredKernelTask>(kernelThreads);

    // PredDivSetup: parallel VELOCITY_FLUX_KERNEL (+ sequential CC_VELOCITY_BC if CC_IBM)
    // Block decomposition: if no Coriolis/patch/CTRL/wind/periodic, use K-block parallel
    bool ccIBM = fds_is_cc_ibm() != 0;
    bool canBlockFlux = fds_velocity_flux_can_block_decompose(1) != 0;
    auto predDivSetupKernelTask = std::make_shared<DivSetupKernelTask>(kernelThreads);

    // WallBC sub-graph: K-block decomposition or mesh-level fallback
    bool canBlockWallBC = fds_wall_bc_can_block_decompose() != 0;
    auto predWallBCSubgraph = canBlockWallBC
        ? buildWallBCBlockSubgraph(nmeshes, blockThreads, numBlocks)
        : buildWallBCSubgraph(nmeshes, kernelThreads);

    // PredWallDiv: parallel PARTICLE_MOMENTUM + DIV_PART_1 kernels
    auto predWallDivKernelTask = std::make_shared<PredWallDivKernelTask>(kernelThreads);

    // PredDivPart2: parallel DIVERGENCE_PART_2_KERNEL
    // Block decomposition: if non-CC_IBM, use K-block parallel
    bool canBlockDivP2 = fds_divergence_part_2_can_block_decompose() != 0 && numBlocks > 1;
    auto predDivP2KernelTask = std::make_shared<DivergencePart2KernelTask>(kernelThreads);

    // VelocityPredictor: block-decomposed kernel (+ CC post-processing collector if CC_IBM)
    // For CC_IBM: skip CFL check in kernel (runs later in collector after CC_PROJECT_VELOCITY)
    auto velPredSubgraph = buildVelocityPredictorBlockSubgraph(
        blockThreads, numBlocks, /*skipCFL=*/ccIBM);
    // Fallback: original mesh-level kernel task for CC_IBM path
    auto velPredKernelTask = std::make_shared<VelocityPredictorKernelTask>(
        kernelThreads, /*skipCFL=*/ccIBM);

    // PredFinal sub-graph (Pattern B, outputs BarrierData)
    auto predFinalSubgraph = buildPredFinalSubgraph(nmeshes, kernelThreads, blockThreads, numBlocks);

    // ChangeTimeStep sub-graph (CFL retry loop)
    auto changeTimeStepSubgraph = buildChangeTimeStepSubgraph(tEnd, nmeshes, kernelThreads);

    // --- Barrier tasks ---

    auto collector1SM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(1)");
    auto meshExchange1 = std::make_shared<MeshExchangeTask>(1, /*ccDensity=*/ccIBM);

    // Merged: HVAC + InitDivIntegrals (eliminates PredInitDivCollector + InitDivTask)
    auto predHvacCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "PredHvacCollector");
    auto hvacInitDivTask = std::make_shared<HvacInitDivTask>(1);

    auto predDivCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "PredDivCollector");
    auto predDivExchangeTask = std::make_shared<DivergenceExchangeTask>(/*corrector=*/false);

    auto predPressureCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "PredPressureCollector");

    auto changeTimeStepCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "ChangeTimeStepCollector");

    auto collector3SM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(3)");
    auto meshExchange3 = std::make_shared<MeshExchangeTask>(3, /*ccDensity=*/false, /*ccEndStep=*/ccIBM);

    // PhaseTransition receives BarrierData directly from PredFinal (no collector needed)
    auto phaseTransTask = std::make_shared<PhaseTransitionTask>();

    // --- Wire the sub-graph ---

    subgraph->inputs(predStep1OrchSM);

    // PredStep1: orchestrator (INSERT_ALL_PARTICLES) -> viscosity -> mass_fd -> DensityPred
    if (canBlockVisc) {
        // Block-decomposed viscosity -> separate mass_fd task
        auto predViscBlockSubgraph = buildComputeViscosityBlockSubgraph(
            nmeshes, blockThreads, numBlocks);
        auto predMassFDKernelTask = std::make_shared<MassFDKernelTask>(kernelThreads);
        subgraph->edges(predStep1OrchSM, predViscBlockSubgraph);
        subgraph->edges(predViscBlockSubgraph, predMassFDKernelTask);
        subgraph->edges(predMassFDKernelTask, densPredKernelTask);
    } else {
        // Mesh-level fallback: combined viscosity + mass_fd
        subgraph->edges(predStep1OrchSM, predStep1KernelTask);
        subgraph->edges(predStep1KernelTask, densPredKernelTask);
    }

    // DensityPred -> MESH_EXCHANGE(1)
    subgraph->edges(densPredKernelTask, collector1SM);
    subgraph->edges(collector1SM, meshExchange1);

    // PredDivSetup: block-decomposed or mesh-level depending on feature flags
    if (canBlockFlux) {
        // Block decomposition: orchestrator(pre-proc + K-decompose) -> parallel blocks -> collector
        // CC_IBM handled internally: CC_VELOCITY_BC + CUTFACE_VELOCITIES in orchestrator,
        // CC_VELOCITY_FLUX in collector
        auto predDivSetupBlockSubgraph = buildVelocityFluxBlockSubgraph(
            nmeshes, blockThreads, numBlocks);
        subgraph->edges(meshExchange1, predDivSetupBlockSubgraph);
        subgraph->edges(predDivSetupBlockSubgraph, predHvacCollectorSM);
    } else if (ccIBM) {
        // CC_IBM with features preventing block decomposition: orchestrator + mesh-level kernel
        auto predDivSetupOrchSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
            std::make_shared<PredDivSetupOrchestrator>(nmeshes), "PredDivSetupOrch");
        subgraph->edges(meshExchange1, predDivSetupOrchSM);
        subgraph->edges(predDivSetupOrchSM, predDivSetupKernelTask);
        subgraph->edges(predDivSetupKernelTask, predHvacCollectorSM);
    } else {
        // Mesh-level fallback (Coriolis, patch velocity, etc.)
        subgraph->edges(meshExchange1, predDivSetupKernelTask);
        subgraph->edges(predDivSetupKernelTask, predHvacCollectorSM);
    }

    // Merged HVAC+InitDiv (was: hvac -> collect -> initDiv)
    subgraph->edges(predHvacCollectorSM, hvacInitDivTask);

    // WallBC sub-graph (three-phase decomposition, reuses corrector pattern)
    subgraph->edges(hvacInitDivTask, predWallBCSubgraph);

    // PredWallDiv: parallel PARTICLE_MOMENTUM + DIV_PART_1 -> DivExchange
    subgraph->edges(predWallBCSubgraph, predWallDivKernelTask);
    subgraph->edges(predWallDivKernelTask, predDivCollectorSM);
    subgraph->edges(predDivCollectorSM, predDivExchangeTask);

    // PredDivPart2 -> Pressure (block-decomposed or mesh-level)
    if (canBlockDivP2) {
        auto predDivP2BlockSubgraph = buildDivergencePart2BlockSubgraph(
            nmeshes, blockThreads, numBlocks);
        subgraph->edges(predDivExchangeTask, predDivP2BlockSubgraph);
        subgraph->edges(predDivP2BlockSubgraph, predPressureCollectorSM);
    } else {
        subgraph->edges(predDivExchangeTask, predDivP2KernelTask);
        subgraph->edges(predDivP2KernelTask, predPressureCollectorSM);
    }

    // Pressure iteration: parallel sub-graph or sequential fallback
    bool useParallelPressure = fds_use_pressure_subgraph() != 0;
    if (useParallelPressure) {
        auto predPressureSubgraph = buildPressureIterationSubgraph(
            tEnd, nmeshes, kernelThreads, /*predictor=*/true, termSignal,
            fds_get_pres_flag());
        subgraph->edges(predPressureCollectorSM, predPressureSubgraph);
        // CC_IBM is always false when useParallelPressure is true
        subgraph->edges(predPressureSubgraph, velPredSubgraph);
    } else {
        auto predPressureTask = std::make_shared<PressureIterationTask>(/*predictor=*/true);
        subgraph->edges(predPressureCollectorSM, predPressureTask);
        if (ccIBM) {
            subgraph->edges(predPressureTask, velPredKernelTask);
        } else {
            subgraph->edges(predPressureTask, velPredSubgraph);
        }
    }

    // VelocityPredictor: block sub-graph or CC_IBM mesh-level path
    if (ccIBM) {
        auto velPredCCSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
            std::make_shared<VelocityPredictorCCCollector>(nmeshes), "VelPredCCCollector");
        subgraph->edges(velPredKernelTask, velPredCCSM);
        subgraph->edges(velPredCCSM, changeTimeStepCollectorSM);
    } else {
        subgraph->edges(velPredSubgraph, changeTimeStepCollectorSM);
    }
    subgraph->edges(changeTimeStepCollectorSM, changeTimeStepSubgraph);
    subgraph->edges(changeTimeStepSubgraph, collector3SM);
    subgraph->edges(collector3SM, meshExchange3);

    // PredFinal (outputs BarrierData) -> PhaseTransition (no collector needed)
    subgraph->edges(meshExchange3, predFinalSubgraph);
    subgraph->edges(predFinalSubgraph, phaseTransTask);

    subgraph->outputs(phaseTransTask);

    return subgraph;
}

#endif // PREDICTOR_SUBGRAPH_H
