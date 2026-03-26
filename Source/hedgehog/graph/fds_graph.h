#ifndef FDS_GRAPH_H
#define FDS_GRAPH_H

#include <hedgehog/hedgehog.h>
#include <service/comm_service.hpp>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/termination_data.h"
#include "../state/timestep_state.h"
#include "../tool/mesh_dependency_graph.h"
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
/// TerminationData is accepted as a second input type and routed to the
/// predictor and corrector sub-graphs, which forward it to their pressure
/// iteration sub-graphs for cycle termination.
///
/// @param nmeshes Number of meshes
/// @param t Initial simulation time
/// @param dt Initial time step
/// @param tEnd End time
/// @param kernelThreads Number of threads for parallel kernel tasks (1 for sequential)
/// @param exchangeThreads Number of threads for parallel flux exchange (default: 1)
/// @return Shared pointer to the constructed graph
inline auto buildFDSGraph(int nmeshes, double t, double dt, double tEnd, size_t kernelThreads,
                          hh::comm::CommService *commService = nullptr,
                          size_t exchangeThreads = 1) {

    using GraphType = hh::Graph<2, MeshData, TerminationData, BarrierData>;
    auto graph = std::make_shared<GraphType>("FDS Hedgehog Graph");

    // --- Build mesh dependency graph (once, shared by predictor & corrector) ---
    auto depGraph = std::make_shared<MeshDependencyGraph>(
        fds_get_lower_mesh_index(), fds_get_lower_mesh_index() + nmeshes - 1);
    // --- Create phase sub-graphs ---
    auto predictorSubgraph = buildPredictorSubgraph(nmeshes, tEnd, kernelThreads, depGraph, commService, exchangeThreads);
    auto correctorSubgraph = buildCorrectorSubgraph(nmeshes, kernelThreads, depGraph, commService, exchangeThreads);

    // --- Create timestep pipeline components ---

    auto timestepDumpSM = std::make_shared<TimestepDumpStateManager>(
        std::make_shared<TimestepDumpState>(nmeshes, tEnd), "TimestepDump");

    auto timestepLoopSM = std::make_shared<TimestepLoopStateManager>(
        std::make_shared<TimestepLoopState>(), "TimestepLoop");
    auto terminationSinkSM = std::make_shared<hh::StateManager<1, BarrierData, BarrierData>>(
        std::make_shared<TerminationSinkState>(), "TerminationSink");

    // --- Wire the graph ---

    // Input: MeshData -> Predictor, TerminationData -> Predictor + Corrector
    graph->inputs(predictorSubgraph);
    graph->input<TerminationData>(correctorSubgraph);

    // Predictor -> Corrector
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
