#ifndef FDS_GRAPH_H
#define FDS_GRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/termination_signal.h"
#include "../state/timestep_state.h"
#include "predictor_subgraph.h"
#include "corrector_subgraph.h"

/// Build the FDS Hedgehog dataflow graph.
///
/// The graph implements the FDS time-stepping loop as a dataflow pipeline:
///   Predictor sub-graph -> Corrector sub-graph -> Timestep loop -> cycle back
///
/// The Corrector sub-graph outputs BarrierData directly (CorrFinal collector
/// already gathers all meshes), so no external TimestepCollector is needed.
///
/// @param nmeshes Number of meshes
/// @param t Initial simulation time
/// @param dt Initial time step
/// @param tEnd End time
/// @param kernelThreads Number of threads for parallel kernel tasks (1 for sequential)
/// @return Shared pointer to the constructed graph
inline auto buildFDSGraph(int nmeshes, double t, double dt, double tEnd, size_t kernelThreads,
                          size_t blockThreads, int numBlocks) {

    using GraphType = hh::Graph<1, MeshData, BarrierData>;
    auto graph = std::make_shared<GraphType>("FDS Hedgehog Graph");

    // --- Shared termination signal for sub-graph cycle termination ---
    auto termSignal = std::make_shared<TerminationSignal>();

    // --- Create phase sub-graphs ---
    auto predictorSubgraph = buildPredictorSubgraph(nmeshes, tEnd, kernelThreads, blockThreads, numBlocks, termSignal);
    auto correctorSubgraph = buildCorrectorSubgraph(nmeshes, tEnd, kernelThreads, blockThreads, numBlocks, termSignal);

    // --- Create timestep pipeline components ---

    auto timestepDumpSM = std::make_shared<TimestepDumpStateManager>(
        std::make_shared<TimestepDumpState>(nmeshes, tEnd), "TimestepDump");

    auto timestepLoopSM = std::make_shared<TimestepLoopStateManager>(
        std::make_shared<TimestepLoopState>(termSignal), "TimestepLoop");
    auto terminationSinkSM = std::make_shared<hh::StateManager<1, BarrierData, BarrierData>>(
        std::make_shared<TerminationSinkState>(), "TerminationSink");

    // --- Wire the graph ---

    // Input -> Predictor -> Corrector
    graph->inputs(predictorSubgraph);
    graph->edges(predictorSubgraph, correctorSubgraph);

    // Corrector -> TimestepDump -> TimestepLoop
    graph->edges(correctorSubgraph, timestepDumpSM);
    graph->edges(timestepDumpSM, timestepLoopSM);

    // Cycle: TimestepLoop -> back to Predictor (MeshData)
    graph->edges(timestepLoopSM, predictorSubgraph);

    // Termination: TimestepLoop -> TerminationSink (BarrierData)
    graph->edges(timestepLoopSM, terminationSinkSM);

    // Graph output: final BarrierData for clean shutdown
    graph->outputs(terminationSinkSM);

    return graph;
}

#endif // FDS_GRAPH_H
