#ifndef FDS_GRAPH_H
#define FDS_GRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/velocity_corrector_data.h"
#include "../task/predictor_tasks.h"
#include "../task/corrector_tasks.h"
#include "../task/barrier_tasks.h"
#include "../task/velocity_corrector_kernel_task.h"
#include "../state/collector_state.h"
#include "../state/mesh_barrier_state.h"
#include "../state/timestep_state.h"
#include "../state/velocity_corrector_state.h"

/// Build the FDS Hedgehog dataflow graph.
///
/// The graph implements the FDS time-stepping loop as a dataflow pipeline:
///   Predictor tasks -> barriers -> Corrector tasks -> barriers -> cycle back
///
/// Barrier pattern: each synchronization point is a CollectorState (pure
/// data-flow: collects N MeshData -> emits 1 BarrierData) followed by a
/// barrier task (computation: receives BarrierData -> emits N MeshData).
///
/// @param nmeshes Number of meshes
/// @param t Initial simulation time
/// @param dt Initial time step
/// @param tEnd End time
/// @param velCorrKernelThreads Number of threads for velocity corrector kernel (1 for sequential)
/// @return Shared pointer to the constructed graph
inline auto buildFDSGraph(int nmeshes, double t, double dt, double tEnd, size_t velCorrKernelThreads) {

    using GraphType = hh::Graph<1, MeshData, MeshData>;
    auto graph = std::make_shared<GraphType>("FDS Hedgehog Graph");

    // --- Create predictor tasks (all sequential) ---
    auto predStep1       = std::make_shared<PredStep1Task>(1);
    auto densityPred     = std::make_shared<DensityPredTask>(1);
    auto predDivSetup    = std::make_shared<PredDivSetupTask>(1);
    auto predWallDiv     = std::make_shared<PredWallDivTask>(1);
    auto divPart2Pred    = std::make_shared<DivPart2PredTask>(1);
    auto velPredictor    = std::make_shared<VelPredictorTask>(1);
    auto predFinal       = std::make_shared<PredFinalTask>(1);

    // --- Create corrector tasks (all sequential) ---
    auto corrStep1       = std::make_shared<CorrStep1Task>(1);
    auto corrDivSetup    = std::make_shared<CorrDivSetupTask>(1);
    auto corrCondens     = std::make_shared<CorrCondensTask>(1);
    auto corrParticle    = std::make_shared<CorrParticleTask>(1);
    auto corrWallBC      = std::make_shared<CorrWallBCTask>(1);
    auto corrRadiation   = std::make_shared<CorrRadiationTask>(1);
    auto corrDivPart1    = std::make_shared<CorrDivPart1Task>(1);
    auto corrDivPart2    = std::make_shared<CorrDivPart2Task>(1);
    // NOTE: corrVelocity replaced by velocity corrector sub-graph (see below)
    auto corrFinal       = std::make_shared<CorrFinalTask>(1);

    // --- Create velocity corrector sub-graph components ---
    // Only the kernel task is parallelized; orchestrator and collector are always sequential
    auto velCorrOrchSM = std::make_shared<hh::StateManager<1, MeshData, VelocityCorrectorWork>>(
        std::make_shared<VelocityCorrectorOrchestrator>(nmeshes), "VelCorrOrch");
    auto velCorrKernelTask = std::make_shared<VelocityCorrectorKernelTask>(velCorrKernelThreads);
    auto velCorrCollectorSM = std::make_shared<hh::StateManager<1, VelocityCorrectorWork, MeshData>>(
        std::make_shared<VelocityCorrectorCollector>(nmeshes), "VelCorrCollector");

    // --- Create barrier collector state managers + barrier tasks ---

    // Predictor: MESH_EXCHANGE(1) after density
    auto collector1SM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(1)");
    auto meshExchange1 = std::make_shared<MeshExchangeTask>(1);

    // Predictor: HVAC barrier
    auto predHvacCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "PredHvacCollector");
    auto predHvacTask = std::make_shared<HvacTask>(1);  // first=1

    // Predictor: INITIALIZE_DIVERGENCE_INTEGRALS
    auto predInitDivCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "PredInitDivCollector");
    auto predInitDivTask = std::make_shared<InitDivIntegralsTask>();

    // Predictor: EXCHANGE_DIVERGENCE_INFO
    auto predDivCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "PredDivCollector");
    auto predDivExchangeTask = std::make_shared<DivergenceExchangeTask>(/*corrector=*/false);

    // Predictor: PRESSURE_ITERATION
    auto predPressureCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "PredPressureCollector");
    auto predPressureTask = std::make_shared<PressureIterationTask>(/*predictor=*/true);

    // Predictor: CHANGE_TIME_STEP_LOOP
    auto changeTimeStepCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "ChangeTimeStepCollector");
    auto changeTimeStepTask = std::make_shared<ChangeTimeStepTask>();

    // Predictor: MESH_EXCHANGE(3) after CFL check
    auto collector3SM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(3)");
    auto meshExchange3 = std::make_shared<MeshExchangeTask>(3);

    // Predictor->Corrector phase transition
    auto phaseTransCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "PhaseTransCollector");
    auto phaseTransTask = std::make_shared<PhaseTransitionTask>();

    // Corrector: MESH_EXCHANGE(4)
    auto collector4SM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(4)");
    auto meshExchange4 = std::make_shared<MeshExchangeTask>(4);

    // Corrector: COMBUSTION barrier
    auto combustionCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "CombustionCollector");
    auto combustionTask = std::make_shared<CombustionTask>();

    // Corrector: HVAC barrier
    auto corrHvacCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "CorrHvacCollector");
    auto corrHvacTask = std::make_shared<HvacTask>(1);  // first=1

    // Corrector: MESH_EXCHANGE(7) particles
    auto collector7SM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(7)");
    auto meshExchange7 = std::make_shared<MeshExchangeTask>(7);

    // Corrector: MESH_EXCHANGE(6) after wall BC
    auto collector6aSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(6a)");
    auto meshExchange6a = std::make_shared<MeshExchangeTask>(6);

    // Corrector: MESH_EXCHANGE(2) after radiation
    auto collector2SM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(2)");
    auto meshExchange2 = std::make_shared<MeshExchangeTask>(2);

    // Corrector: INITIALIZE_DIVERGENCE_INTEGRALS
    auto corrInitDivCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "CorrInitDivCollector");
    auto corrInitDivTask = std::make_shared<InitDivIntegralsTask>();

    // Corrector: EXCHANGE_DIVERGENCE_INFO + RTE
    auto corrDivCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "CorrDivCollector");
    auto corrDivExchangeTask = std::make_shared<DivergenceExchangeTask>(/*corrector=*/true);

    // Corrector: PRESSURE_ITERATION
    auto corrPressureCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "CorrPressureCollector");
    auto corrPressureTask = std::make_shared<PressureIterationTask>(/*predictor=*/false);

    // Corrector: MESH_EXCHANGE(6) after velocity
    auto collector6bSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "Collector(6b)");
    auto meshExchange6b = std::make_shared<MeshExchangeTask>(6);

    // Passthrough barriers (pure data flow, no computation — kept as direct states)
    auto predStep1BarrierSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<PassthroughBarrierState>(nmeshes), "PredStep1Barrier");
    auto corrCondensBarrierSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<PassthroughBarrierState>(nmeshes), "CorrCondensBarrier");

    // Timestep: collector -> computation task -> loop state (cycle management)
    auto timestepCollectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CollectorState>(nmeshes), "TimestepCollector");
    auto timestepTask = std::make_shared<TimestepTask>(tEnd);
    auto timestepLoopState = std::make_shared<TimestepLoopState>();
    auto timestepLoopSM = std::make_shared<TimestepLoopStateManager>(timestepLoopState);

    // --- Wire the graph ---

    // Graph input goes to predStep1
    graph->inputs(predStep1);

    // Predictor pipeline
    graph->edges(predStep1, predStep1BarrierSM);             // Passthrough barrier
    graph->edges(predStep1BarrierSM, densityPred);
    graph->edges(densityPred, collector1SM);                  // Collect for MESH_EXCHANGE(1)
    graph->edges(collector1SM, meshExchange1);                // Do MESH_EXCHANGE(1)
    graph->edges(meshExchange1, predDivSetup);
    graph->edges(predDivSetup, predHvacCollectorSM);          // Collect for HVAC
    graph->edges(predHvacCollectorSM, predHvacTask);          // Do HVAC_CALC
    graph->edges(predHvacTask, predInitDivCollectorSM);       // Collect for INIT_DIV
    graph->edges(predInitDivCollectorSM, predInitDivTask);    // Do INIT_DIV_INTEGRALS
    graph->edges(predInitDivTask, predWallDiv);               // wall_bc + particle_momentum + div_part_1
    graph->edges(predWallDiv, predDivCollectorSM);            // Collect for DIV_EXCHANGE
    graph->edges(predDivCollectorSM, predDivExchangeTask);    // Do EXCHANGE_DIV_INFO
    graph->edges(predDivExchangeTask, divPart2Pred);
    graph->edges(divPart2Pred, predPressureCollectorSM);      // Collect for PRESSURE
    graph->edges(predPressureCollectorSM, predPressureTask);  // Do PRESSURE_ITERATION
    graph->edges(predPressureTask, velPredictor);
    graph->edges(velPredictor, changeTimeStepCollectorSM);    // Collect for CFL check
    graph->edges(changeTimeStepCollectorSM, changeTimeStepTask); // Do CHANGE_TIME_STEP_LOOP
    graph->edges(changeTimeStepTask, collector3SM);            // Collect for MESH_EXCHANGE(3)
    graph->edges(collector3SM, meshExchange3);                 // Do MESH_EXCHANGE(3)
    graph->edges(meshExchange3, predFinal);
    graph->edges(predFinal, phaseTransCollectorSM);           // Collect for phase transition
    graph->edges(phaseTransCollectorSM, phaseTransTask);      // Do phase transition
    graph->edges(phaseTransTask, corrStep1);                  // -> corrector

    // Corrector pipeline
    graph->edges(corrStep1, collector4SM);                    // Collect for MESH_EXCHANGE(4)
    graph->edges(collector4SM, meshExchange4);                // Do MESH_EXCHANGE(4)
    graph->edges(meshExchange4, corrDivSetup);
    graph->edges(corrDivSetup, combustionCollectorSM);        // Collect for COMBUSTION
    graph->edges(combustionCollectorSM, combustionTask);      // Do COMBUSTION
    graph->edges(combustionTask, corrHvacCollectorSM);        // Collect for HVAC
    graph->edges(corrHvacCollectorSM, corrHvacTask);          // Do HVAC_CALC
    graph->edges(corrHvacTask, corrCondens);
    graph->edges(corrCondens, corrCondensBarrierSM);          // Passthrough barrier
    graph->edges(corrCondensBarrierSM, corrParticle);
    graph->edges(corrParticle, collector7SM);                 // Collect for MESH_EXCHANGE(7)
    graph->edges(collector7SM, meshExchange7);                // Do MESH_EXCHANGE(7)
    graph->edges(meshExchange7, corrWallBC);
    graph->edges(corrWallBC, collector6aSM);                  // Collect for MESH_EXCHANGE(6)
    graph->edges(collector6aSM, meshExchange6a);              // Do MESH_EXCHANGE(6)
    graph->edges(meshExchange6a, corrRadiation);
    graph->edges(corrRadiation, collector2SM);                // Collect for MESH_EXCHANGE(2)
    graph->edges(collector2SM, meshExchange2);                // Do MESH_EXCHANGE(2)
    graph->edges(meshExchange2, corrInitDivCollectorSM);      // Collect for INIT_DIV
    graph->edges(corrInitDivCollectorSM, corrInitDivTask);    // Do INIT_DIV_INTEGRALS
    graph->edges(corrInitDivTask, corrDivPart1);
    graph->edges(corrDivPart1, corrDivCollectorSM);           // Collect for DIV_EXCHANGE + RTE
    graph->edges(corrDivCollectorSM, corrDivExchangeTask);    // Do EXCHANGE_DIV_INFO + RTE
    graph->edges(corrDivExchangeTask, corrDivPart2);
    graph->edges(corrDivPart2, corrPressureCollectorSM);      // Collect for PRESSURE
    graph->edges(corrPressureCollectorSM, corrPressureTask);  // Do PRESSURE_ITERATION
    // Velocity corrector sub-graph (parallel multi-mesh execution)
    graph->edges(corrPressureTask, velCorrOrchSM);            // Pressure → VelCorrOrchestrator
    graph->edges(velCorrOrchSM, velCorrKernelTask);           // Orchestrator → Kernel (parallel)
    graph->edges(velCorrKernelTask, velCorrCollectorSM);      // Kernel → Collector
    graph->edges(velCorrCollectorSM, collector6bSM);          // Collector → MESH_EXCHANGE(6)
    graph->edges(collector6bSM, meshExchange6b);              // Do MESH_EXCHANGE(6)
    graph->edges(meshExchange6b, corrFinal);

    // End of time step: corrFinal -> collector -> timestep task -> loop state -> cycle
    graph->edges(corrFinal, timestepCollectorSM);
    graph->edges(timestepCollectorSM, timestepTask);
    graph->edges(timestepTask, timestepLoopSM);

    // Cycle: timestep loop state -> back to predictor
    graph->edges(timestepLoopSM, predStep1);

    // Graph output (for termination detection)
    graph->outputs(timestepLoopSM);

    return graph;
}

#endif // FDS_GRAPH_H
