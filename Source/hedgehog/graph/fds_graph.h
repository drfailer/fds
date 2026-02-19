#ifndef FDS_GRAPH_H
#define FDS_GRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../task/predictor_tasks.h"
#include "../task/corrector_tasks.h"
#include "../state/mesh_barrier_state.h"
#include "../state/divergence_barrier_state.h"
#include "../state/pressure_barrier_state.h"
#include "../state/phase_transition_state.h"
#include "../state/timestep_state.h"

/// Build the FDS Hedgehog dataflow graph.
///
/// The graph implements the FDS time-stepping loop as a dataflow pipeline:
///   Predictor tasks -> barriers -> Corrector tasks -> barriers -> cycle back
///
/// Phase 1: numThreads=1 for all tasks (sequential, for correctness verification)
/// Phase 2: numThreads=nmeshes for per-mesh tasks (parallel mesh processing)
///
/// @param nmeshes Number of meshes
/// @param t Initial simulation time
/// @param dt Initial time step
/// @param tEnd End time
/// @param numThreads Number of threads per task (1 for Phase 1)
/// @return Shared pointer to the constructed graph
inline auto buildFDSGraph(int nmeshes, double t, double dt, double tEnd, size_t numThreads) {

    using GraphType = hh::Graph<1, MeshData, MeshData>;
    auto graph = std::make_shared<GraphType>("FDS Hedgehog Graph");

    // --- Create predictor tasks ---
    auto predStep1       = std::make_shared<PredStep1Task>(numThreads);
    auto densityPred     = std::make_shared<DensityPredTask>(numThreads);
    auto predDivSetup    = std::make_shared<PredDivSetupTask>(numThreads);
    auto predWallDiv     = std::make_shared<PredWallDivTask>(numThreads);
    auto divPart2Pred    = std::make_shared<DivPart2PredTask>(numThreads);
    auto velPredictor    = std::make_shared<VelPredictorTask>(numThreads);
    auto predFinal       = std::make_shared<PredFinalTask>(numThreads);

    // --- Create corrector tasks ---
    auto corrStep1       = std::make_shared<CorrStep1Task>(numThreads);
    auto corrDivSetup    = std::make_shared<CorrDivSetupTask>(numThreads);
    auto corrCondens     = std::make_shared<CorrCondensTask>(numThreads);
    auto corrParticle    = std::make_shared<CorrParticleTask>(numThreads);
    auto corrWallBC      = std::make_shared<CorrWallBCTask>(numThreads);
    auto corrRadiation   = std::make_shared<CorrRadiationTask>(numThreads);
    auto corrDivPart1    = std::make_shared<CorrDivPart1Task>(numThreads);
    auto corrDivPart2    = std::make_shared<CorrDivPart2Task>(numThreads);
    auto corrVelocity    = std::make_shared<CorrVelocityTask>(numThreads);
    auto corrFinal       = std::make_shared<CorrFinalTask>(numThreads);

    // --- Create barrier state managers ---

    // Predictor barriers
    auto barrier1SM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<MeshBarrierState>(nmeshes, 1), "Barrier(1)");

    // Initialize divergence integrals BEFORE divergence_part_1 (predictor)
    auto predInitDivSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<InitDivIntegralsBarrier>(nmeshes), "PredInitDiv");

    // Exchange divergence info AFTER divergence_part_1 (predictor)
    auto predDivBarrierSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<DivergenceBarrierState>(nmeshes, /*corrector=*/false), "PredDivBarrier");

    auto predPressureBarrierSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<PressureBarrierState>(nmeshes), "PredPressureBarrier");

    auto barrier3SM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<MeshBarrierState>(nmeshes, 3), "Barrier(3)");

    // Predictor->Corrector transition
    auto phaseTransSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<PhaseTransitionState>(nmeshes), "PhaseTransition");

    // Corrector barriers
    auto barrier4SM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<MeshBarrierState>(nmeshes, 4), "Barrier(4)");

    auto corrCombustionBarrierSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<MeshBarrierState>(nmeshes, 0), "CombustionBarrier");
    // Note: combustion is load-balanced across all meshes, so it acts as a barrier

    auto barrier7SM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<MeshBarrierState>(nmeshes, 7), "Barrier(7)");

    auto barrier6aSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<MeshBarrierState>(nmeshes, 6), "Barrier(6a)");

    auto barrier2SM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<MeshBarrierState>(nmeshes, 2), "Barrier(2)");

    // Initialize divergence integrals BEFORE divergence_part_1 (corrector)
    auto corrInitDivSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<InitDivIntegralsBarrier>(nmeshes), "CorrInitDiv");

    // Exchange divergence info AFTER divergence_part_1 (corrector)
    // Also calls RTE source correction in corrector phase
    auto corrDivBarrierSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<DivergenceBarrierState>(nmeshes, /*corrector=*/true), "CorrDivBarrier");

    auto corrPressureBarrierSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<PressureBarrierState>(nmeshes), "CorrPressureBarrier");

    auto barrier6bSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<MeshBarrierState>(nmeshes, 6), "Barrier(6b)");

    // Timestep loop state (with cycle detection)
    auto timestepState = std::make_shared<TimestepState>(nmeshes, tEnd);
    auto timestepSM = std::make_shared<TimestepStateManager>(timestepState);

    // --- Wire the graph ---

    // Graph input goes to predStep1
    graph->inputs(predStep1);

    // Predictor pipeline
    graph->edges(predStep1, densityPred);
    graph->edges(densityPred, barrier1SM);              // Barrier: MESH_EXCHANGE(1)
    graph->edges(barrier1SM, predDivSetup);
    graph->edges(predDivSetup, predInitDivSM);          // Barrier: INITIALIZE_DIVERGENCE_INTEGRALS
    graph->edges(predInitDivSM, predWallDiv);           // wall_bc + particle_momentum + divergence_part_1
    graph->edges(predWallDiv, predDivBarrierSM);        // Barrier: EXCHANGE_DIVERGENCE_INFO
    graph->edges(predDivBarrierSM, divPart2Pred);
    graph->edges(divPart2Pred, predPressureBarrierSM);  // Barrier: PRESSURE_ITERATION_SCHEME
    graph->edges(predPressureBarrierSM, velPredictor);
    graph->edges(velPredictor, barrier3SM);              // Barrier: MESH_EXCHANGE(3)
    graph->edges(barrier3SM, predFinal);
    graph->edges(predFinal, phaseTransSM);              // Barrier: Phase transition

    // Corrector pipeline
    graph->edges(phaseTransSM, corrStep1);
    graph->edges(corrStep1, barrier4SM);                // Barrier: MESH_EXCHANGE(4)
    graph->edges(barrier4SM, corrDivSetup);
    graph->edges(corrDivSetup, corrCombustionBarrierSM); // Barrier: combustion
    graph->edges(corrCombustionBarrierSM, corrCondens);
    graph->edges(corrCondens, corrParticle);
    graph->edges(corrParticle, barrier7SM);              // Barrier: MESH_EXCHANGE(7) particles
    graph->edges(barrier7SM, corrWallBC);
    graph->edges(corrWallBC, barrier6aSM);               // Barrier: MESH_EXCHANGE(6) back wall
    graph->edges(barrier6aSM, corrRadiation);
    graph->edges(corrRadiation, barrier2SM);              // Barrier: MESH_EXCHANGE(2) radiation
    graph->edges(barrier2SM, corrInitDivSM);             // Barrier: INITIALIZE_DIVERGENCE_INTEGRALS
    graph->edges(corrInitDivSM, corrDivPart1);           // combustion_bc + divergence_part_1
    graph->edges(corrDivPart1, corrDivBarrierSM);        // Barrier: EXCHANGE_DIVERGENCE_INFO + RTE
    graph->edges(corrDivBarrierSM, corrDivPart2);
    graph->edges(corrDivPart2, corrPressureBarrierSM);   // Barrier: PRESSURE_ITERATION_SCHEME
    graph->edges(corrPressureBarrierSM, corrVelocity);
    graph->edges(corrVelocity, barrier6bSM);             // Barrier: MESH_EXCHANGE(6)
    graph->edges(barrier6bSM, corrFinal);
    graph->edges(corrFinal, timestepSM);                 // End of time step

    // Cycle: timestep state -> back to predictor
    graph->edges(timestepSM, predStep1);

    // Graph output (for termination detection)
    graph->outputs(timestepSM);

    return graph;
}

#endif // FDS_GRAPH_H
