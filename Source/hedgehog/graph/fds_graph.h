#ifndef FDS_GRAPH_H
#define FDS_GRAPH_H

#include <hedgehog/hedgehog.h>
#include <service/comm_service.hpp>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/termination_data.h"
#include "../state/timestep_state.h"
#include "../task/timestep_tasks.h"
#include "../tool/mesh_dependency_graph.h"
#include "../tool/thread_budget.h"
#include "predictor_subgraph.h"
#include "corrector_subgraph.h"

/// Build the FDS Hedgehog dataflow graph.
///
/// The graph implements the FDS time-stepping loop as a dataflow pipeline:
///   Predictor -> Corrector -> Dump fork -> TimestepState -> cycle back
///
/// The dump phase uses a fork-join pattern:
///   PreDumpScatter forks into two parallel branches:
///     - DumpGlobalTask: global computation + global file I/O
///     - DumpMeshOutputsTask: per-mesh file I/O (skipped on non-dump timesteps)
///   TimestepState joins both branches, runs STOP_CHECK, then either cycles
///   MeshData back to the predictor or emits BarrierData for termination.
///
/// @param nmeshes Number of meshes
/// @param t Initial simulation time
/// @param dt Initial time step
/// @param tEnd End time
/// @param budget Thread budget computed from hardware capabilities
/// @return Shared pointer to the constructed graph
inline auto buildFDSGraph(int nmeshes, double t, double dt, double tEnd,
                          const ThreadBudget &budget,
                          hh::comm::CommService *commService = nullptr) {

    using GraphType = hh::Graph<2, MeshData, TerminationData, BarrierData>;
    auto graph = std::make_shared<GraphType>("FDS Hedgehog Graph");

    // --- Build mesh dependency graph (once, shared by predictor & corrector) ---
    auto depGraph = std::make_shared<MeshDependencyGraph>(
        fds_get_lower_mesh_index(), fds_get_lower_mesh_index() + nmeshes - 1);
    // --- Create phase sub-graphs ---
    auto predictorSubgraph = buildPredictorSubgraph(nmeshes, budget, depGraph, commService);
    auto correctorSubgraph = buildCorrectorSubgraph(nmeshes, budget, depGraph, commService);

    // --- Create dump fork-join + timestep loop ---

    // Shared icyc counter between DumpGlobal (SET_DIAGNOSTICS) and TimestepState (SET_ICYC)
    auto icyc = std::make_shared<int>(1);

    // Scatter: BarrierData -> fork into MeshData (per-mesh) + BarrierData (global)
    // Checks dump schedule; skips MeshData emission on non-dump timesteps.
    auto preDumpScatter = std::make_shared<PreDumpScatterTask>();

    // Fork branch 1: global computation + global file I/O (1 thread)
    auto dumpGlobalTask = std::make_shared<DumpGlobalTask>(icyc);

    // Fork branch 2: parallel per-mesh dump I/O (N threads, idle on non-dump steps)
    auto dumpMeshTask = std::make_shared<DumpMeshOutputsTask>(
        static_cast<size_t>(nmeshes));

    // Merged join + timestep loop: collects dump results, STOP_CHECK,
    // then cycles MeshData back or emits BarrierData for termination.
    auto timestepSM = std::make_shared<TimestepStateManager>(
        std::make_shared<TimestepState>(nmeshes, tEnd, icyc), "Timestep");

    // --- Wire the graph ---

    // Input: MeshData -> Predictor, TerminationData -> Predictor + Corrector
    graph->inputs(predictorSubgraph);
    graph->input<TerminationData>(correctorSubgraph);

    // Predictor -> Corrector
    graph->edges(predictorSubgraph, correctorSubgraph);

    // Corrector -> Scatter (BarrierData)
    graph->edges(correctorSubgraph, preDumpScatter);

    // Fork: Scatter -> DumpGlobal (BarrierData) + DumpMesh (MeshData)
    graph->edges(preDumpScatter, dumpGlobalTask);
    graph->edges(preDumpScatter, dumpMeshTask);

    // Join: DumpGlobal (BarrierData) + DumpMesh (MeshData) -> TimestepState
    graph->edges(dumpGlobalTask, timestepSM);
    graph->edges(dumpMeshTask, timestepSM);

    // Cycle: TimestepState -> back to Predictor (MeshData)
    graph->edges(timestepSM, predictorSubgraph);

    // Graph output: TimestepState emits BarrierData when simulation is done
    graph->outputs(timestepSM);

    return graph;
}

#endif // FDS_GRAPH_H
