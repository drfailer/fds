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
///   Predictor -> Corrector -> Dump fork -> Timestep loop -> cycle back
///
/// The dump phase uses a fork-join pattern:
///   PreDumpScatter forks into two parallel branches:
///     - DumpGlobalTask: global computation (SET_DIAGNOSTICS, EXCHANGE_GLOBAL_OUTPUTS,
///       UPDATE_CONTROLS) + global file I/O (HRR.csv, mass.csv, devc.csv, .smv)
///     - DumpMeshOutputsTask: per-mesh file I/O (SLCF .sf, BNDF .bf, PRT5 .prt5, etc.)
///   PostDump joins both branches, then runs STOP_CHECK + termination decision.
///
/// Global files and per-mesh files are completely independent, so the two
/// branches run in parallel with no contention.
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

    // --- Create dump fork-join components ---

    // Shared icyc counter between DumpGlobal (SET_DIAGNOSTICS) and PostDump (SET_ICYC)
    auto icyc = std::make_shared<int>(1);

    // Scatter: BarrierData -> fork into MeshData (per-mesh) + BarrierData (global)
    auto preDumpScatter = std::make_shared<PreDumpScatterTask>();

    // Fork branch 1: global computation + global file I/O (1 thread)
    auto dumpGlobalTask = std::make_shared<DumpGlobalTask>(icyc);

    // Fork branch 2: parallel per-mesh dump I/O (N threads)
    auto dumpMeshTask = std::make_shared<DumpMeshOutputsTask>(
        static_cast<size_t>(nmeshes));

    // Join: collect N MeshData + 1 BarrierData, then STOP_CHECK + termination
    auto postDumpSM = std::make_shared<PostDumpStateManager>(
        std::make_shared<PostDumpState>(nmeshes, tEnd, icyc), "PostDump");

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

    // Corrector -> Scatter (BarrierData)
    graph->edges(correctorSubgraph, preDumpScatter);

    // Fork: Scatter -> DumpGlobal (BarrierData) + DumpMesh (MeshData)
    graph->edges(preDumpScatter, dumpGlobalTask);
    graph->edges(preDumpScatter, dumpMeshTask);

    // Join: DumpGlobal (BarrierData) + DumpMesh (MeshData) -> PostDump
    graph->edges(dumpGlobalTask, postDumpSM);
    graph->edges(dumpMeshTask, postDumpSM);

    // PostDump -> TimestepLoop (BarrierData)
    graph->edges(postDumpSM, timestepLoopSM);

    // Cycle: TimestepLoop -> back to Predictor (MeshData)
    graph->edges(timestepLoopSM, predictorSubgraph);

    // Termination: TimestepLoop -> TerminationSink (BarrierData)
    graph->edges(timestepLoopSM, terminationSinkSM);

    // Graph output: final BarrierData for clean shutdown
    graph->outputs(terminationSinkSM);

    return graph;
}

#endif // FDS_GRAPH_H
