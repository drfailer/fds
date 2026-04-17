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
#include "pressure_iteration_subgraph.h"

/// Build the FDS Hedgehog dataflow graph.
///
/// The graph implements the FDS time-stepping loop as a dataflow pipeline:
///   Predictor -> Corrector -> Dump fork -> TimestepState -> cycle back
///
/// The dump phase uses a fork-join pattern:
///   CorrFinalDumpTask forks into two parallel branches:
///     - DumpGlobalTask (BarrierData): global computation + global file I/O
///     - DumpMeshOutputsTask (MeshData<>): per-mesh file I/O (skipped on non-dump timesteps)
///   TimestepState joins both branches, runs STOP_CHECK, then either cycles
///   MeshData<> back to the predictor or emits BarrierData for termination.
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

    using GraphType = hh::Graph<2, MeshData<MeshState::Init>, TerminationData, BarrierData>;
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

    // Fork branch 1: global computation + global file I/O (1 thread)
    auto dumpGlobalTask = std::make_shared<DumpGlobalTask>(icyc);

    // Fork branch 2: parallel per-mesh dump I/O (N threads, idle on non-dump steps)
    auto dumpMeshTask = std::make_shared<DumpMeshOutputsTask>(
        static_cast<size_t>(nmeshes));

    // Merged join + timestep loop: collects dump results, STOP_CHECK,
    // then cycles MeshData<> back or emits BarrierData for termination.
    auto timestepSM = std::make_shared<TimestepStateManager>(
        std::make_shared<TimestepState>(nmeshes, tEnd, icyc), "Timestep");

    bool useParallelPressure = fds_use_pressure_subgraph() != 0;

    // --- Wire the graph ---

    // Input: MeshData<Init> -> TimestepState (INSERT_PARTICLES, then emit to Predictor)
    graph->input<MeshData<MeshState::Init>>(timestepSM);
    // Input: TerminationData -> Predictor + Corrector (for cycle termination)
    graph->input<TerminationData>(predictorSubgraph);
    graph->input<TerminationData>(correctorSubgraph);

    // Predictor -> Corrector (MeshData<> only; PredPressure routes elsewhere)
    graph->edges(predictorSubgraph, correctorSubgraph);

    // Fork: Corrector -> DumpGlobal (BarrierData) + DumpMesh (MeshData<>)
    // CorrFinalDumpTask checks dump schedule and emits MeshData<> only when needed.
    graph->edges(correctorSubgraph, dumpGlobalTask);
    graph->edges(correctorSubgraph, dumpMeshTask);

    // Join: DumpGlobal (BarrierData) + DumpMesh (MeshData<>) -> TimestepState
    graph->edges(dumpGlobalTask, timestepSM);
    graph->edges(dumpMeshTask, timestepSM);

    // Cycle: TimestepState -> back to Predictor (MeshData<>)
    graph->edges(timestepSM, predictorSubgraph);

    // Graph output: TimestepState emits BarrierData when simulation is done
    graph->outputs(timestepSM);

    // --- Shared pressure subgraph (when enabled) ---
    if (useParallelPressure) {
        auto pressureSubgraph = buildPressureIterationSubgraph(
            nmeshes, budget, depGraph, commService, fds_get_pres_flag());

        graph->input<TerminationData>(pressureSubgraph);

        // Predictor <-> PressureSubgraph (PredPressure)
        graph->edges(predictorSubgraph, pressureSubgraph);
        graph->edges(pressureSubgraph, predictorSubgraph);

        // Corrector <-> PressureSubgraph (CorrPressure)
        graph->edges(correctorSubgraph, pressureSubgraph);
        graph->edges(pressureSubgraph, correctorSubgraph);
    }

    return graph;
}

#endif // FDS_GRAPH_H
