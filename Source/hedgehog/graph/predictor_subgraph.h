#ifndef PREDICTOR_SUBGRAPH_H
#define PREDICTOR_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/velocity_predictor_data.h"
#include "../state/collector_state.h"
#include "../state/pred_step1_state.h"
#include "../state/div_setup_state.h"
#include "../state/velocity_predictor_state.h"
#include "../task/barrier_tasks.h"
#include "../task/pred_step1_kernel_task.h"
#include "../task/density_pred_kernel_task.h"
#include "../task/div_setup_kernel_task.h"
#include "../task/pred_wall_div_kernel_task.h"
#include "../task/divergence_part2_kernel_task.h"
#include "../task/velocity_predictor_kernel_task.h"
#include "change_timestep_subgraph.h"
#include "velocity_bc_subgraph.h"
#include "wallbc_subgraph.h"

/// Build the Predictor sub-graph.
///
/// Implements the full predictor phase of the FDS time-stepping loop:
///   PredStep1 -> DensityPred -> MESH_EXCHANGE(1) -> PredDivSetup -> HVAC ->
///   InitDivIntegrals -> WallBC -> PredWallDiv -> DivergenceExchange -> PredDivPart2 ->
///   PressureIteration -> VelocityPredictor -> ChangeTimeStep ->
///   MESH_EXCHANGE(3) -> PredFinal -> PhaseTransition
///
/// @param nmeshes Number of meshes
/// @param tEnd Simulation end time (passed to ChangeTimeStep sub-graph)
/// @param kernelThreads Number of threads for parallel kernel tasks
/// @return Shared pointer to the constructed sub-graph
inline auto buildPredictorSubgraph(int nmeshes, double tEnd, size_t kernelThreads) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData, MeshData>>("Predictor");

    // --- Kernel sub-graph components ---

    // PredStep1: sequential INSERT_ALL_PARTICLES -> parallel kernels
    auto predStep1OrchSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<PredStep1Orchestrator>(nmeshes), "PredStep1Orch");
    auto predStep1KernelTask = std::make_shared<PredStep1KernelTask>(kernelThreads);

    // DensityPred: parallel DENSITY_KERNEL
    auto densPredKernelTask = std::make_shared<DensityPredKernelTask>(kernelThreads);

    // PredDivSetup: sequential CC_VELOCITY_BC -> parallel VELOCITY_FLUX_KERNEL
    auto predDivSetupOrchSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<PredDivSetupOrchestrator>(nmeshes), "PredDivSetupOrch");
    auto predDivSetupKernelTask = std::make_shared<DivSetupKernelTask>(kernelThreads);

    // WallBC sub-graph (three-phase: preprocessing -> parallel kernel -> finalize)
    auto predWallBCSubgraph = buildWallBCSubgraph(nmeshes, kernelThreads);

    // PredWallDiv: parallel PARTICLE_MOMENTUM + DIV_PART_1 kernels
    auto predWallDivKernelTask = std::make_shared<PredWallDivKernelTask>(kernelThreads);

    // PredDivPart2: parallel DIVERGENCE_PART_2_KERNEL
    auto predDivP2KernelTask = std::make_shared<DivergencePart2KernelTask>(kernelThreads);

    // VelocityPredictor: parallel kernel + CC_PROJECT_VELOCITY
    auto velPredOrchSM = std::make_shared<hh::StateManager<1, MeshData, VelocityPredictorWork>>(
        std::make_shared<VelocityPredictorOrchestrator>(nmeshes), "VelPredOrch");
    auto velPredKernelTask = std::make_shared<VelocityPredictorKernelTask>(kernelThreads);
    auto velPredCollectorSM = std::make_shared<hh::StateManager<1, VelocityPredictorWork, MeshData>>(
        std::make_shared<VelocityPredictorCollector>(nmeshes), "VelPredCollector");

    // PredFinal sub-graph (Pattern B)
    auto predFinalSubgraph = buildPredFinalSubgraph(nmeshes, kernelThreads);

    // ChangeTimeStep sub-graph (CFL retry loop)
    auto changeTimeStepSubgraph = buildChangeTimeStepSubgraph(tEnd);

    // --- Barrier tasks ---

    auto collector1SM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(1)");
    auto meshExchange1 = std::make_shared<MeshExchangeTask>(1);

    auto predHvacCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "PredHvacCollector");
    auto predHvacTask = std::make_shared<HvacTask>(1);

    auto predInitDivCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "PredInitDivCollector");
    auto predInitDivTask = std::make_shared<InitDivIntegralsTask>();

    auto predDivCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "PredDivCollector");
    auto predDivExchangeTask = std::make_shared<DivergenceExchangeTask>(/*corrector=*/false);

    auto predPressureCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "PredPressureCollector");
    auto predPressureTask = std::make_shared<PressureIterationTask>(/*predictor=*/true);

    auto changeTimeStepCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "ChangeTimeStepCollector");

    auto collector3SM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(3)");
    auto meshExchange3 = std::make_shared<MeshExchangeTask>(3);

    auto phaseTransCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "PhaseTransCollector");
    auto phaseTransTask = std::make_shared<PhaseTransitionTask>();

    // --- Wire the sub-graph ---

    subgraph->inputs(predStep1OrchSM);

    // PredStep1: orchestrator (INSERT_ALL_PARTICLES) -> parallel kernel -> DensityPred
    subgraph->edges(predStep1OrchSM, predStep1KernelTask);
    subgraph->edges(predStep1KernelTask, densPredKernelTask);

    // DensityPred -> MESH_EXCHANGE(1)
    subgraph->edges(densPredKernelTask, collector1SM);
    subgraph->edges(collector1SM, meshExchange1);

    // PredDivSetup: orchestrator (CC_VELOCITY_BC) -> parallel kernel -> HVAC
    subgraph->edges(meshExchange1, predDivSetupOrchSM);
    subgraph->edges(predDivSetupOrchSM, predDivSetupKernelTask);
    subgraph->edges(predDivSetupKernelTask, predHvacCollectorSM);
    subgraph->edges(predHvacCollectorSM, predHvacTask);
    subgraph->edges(predHvacTask, predInitDivCollectorSM);
    subgraph->edges(predInitDivCollectorSM, predInitDivTask);

    // WallBC sub-graph (three-phase decomposition, reuses corrector pattern)
    subgraph->edges(predInitDivTask, predWallBCSubgraph);

    // PredWallDiv: parallel PARTICLE_MOMENTUM + DIV_PART_1 -> DivExchange
    subgraph->edges(predWallBCSubgraph, predWallDivKernelTask);
    subgraph->edges(predWallDivKernelTask, predDivCollectorSM);
    subgraph->edges(predDivCollectorSM, predDivExchangeTask);

    // PredDivPart2 -> Pressure
    subgraph->edges(predDivExchangeTask, predDivP2KernelTask);
    subgraph->edges(predDivP2KernelTask, predPressureCollectorSM);
    subgraph->edges(predPressureCollectorSM, predPressureTask);

    // VelocityPredictor (keeps orch/collector for CC_PROJECT_VELOCITY)
    subgraph->edges(predPressureTask, velPredOrchSM);
    subgraph->edges(velPredOrchSM, velPredKernelTask);
    subgraph->edges(velPredKernelTask, velPredCollectorSM);
    subgraph->edges(velPredCollectorSM, changeTimeStepCollectorSM);
    subgraph->edges(changeTimeStepCollectorSM, changeTimeStepSubgraph);
    subgraph->edges(changeTimeStepSubgraph, collector3SM);
    subgraph->edges(collector3SM, meshExchange3);

    // PredFinal sub-graph + phase transition
    subgraph->edges(meshExchange3, predFinalSubgraph);
    subgraph->edges(predFinalSubgraph, phaseTransCollectorSM);
    subgraph->edges(phaseTransCollectorSM, phaseTransTask);

    subgraph->outputs(phaseTransTask);

    return subgraph;
}

#endif // PREDICTOR_SUBGRAPH_H
