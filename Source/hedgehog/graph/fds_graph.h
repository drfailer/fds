#ifndef FDS_GRAPH_H
#define FDS_GRAPH_H

#include <hedgehog/hedgehog.h>
#include <hedgehog_comm.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/termination_data.h"
#include "../task/timestep_tasks.h"
#include "../tool/thread_budget.h"
#include "../exchange/mesh_dependency_graph.h"
#include "../exchange/exchange_graph.h"
#include "../exchange/exchange_dispatch.h"
#include "predictor_subgraph.h"
#include "corrector_subgraph.h"
#include "pressure_iteration_subgraph.h"

/// Build the FDS Hedgehog dataflow graph.
///
/// The graph implements the FDS time-stepping loop as a dataflow pipeline:
///   Predictor -> Corrector -> TimestepTask -> cycle back
///
/// TimestepTask handles end-of-corrector work (RTE, reduce HRR/mass),
/// dump I/O (global via AsyncWorker || per-mesh via ThreadPool), STOP_CHECK,
/// and DT adjustment. Emits BarrierData for termination.
///
/// @param nmeshes Number of meshes
/// @param t Initial simulation time
/// @param dt Initial time step
/// @param tEnd End time
/// @param budget Thread budget computed from hardware capabilities
/// @return Shared pointer to the constructed graph
inline auto buildFDSGraph(int nmeshes, double t, double dt, double tEnd,
                          const ThreadBudget &budget,
                          std::shared_ptr<MeshDependencyGraph> depGraph = nullptr,
                          hh::comm::CommService *commService = nullptr) {

    using GraphType = hh::Graph<2, MeshData<MeshState::Init>, TerminationData, BarrierData>;
    auto graph = std::make_shared<GraphType>("FDS Hedgehog Graph");

    // --- Create phase sub-graphs ---
    auto predictorSubgraph = buildPredictorSubgraph(nmeshes, budget);
    auto correctorSubgraph = buildCorrectorSubgraph(nmeshes, budget);

    // --- Create timestep task ---

    auto icyc = std::make_shared<int>(1);

    // Merged timestep: RTE, reduce, dump (pool+async), stop_check, cycle.
    auto timestepTask = std::make_shared<TimestepTask>(
        nmeshes, static_cast<size_t>(nmeshes), tEnd, icyc);

    bool useParallelPressure = fds_use_pressure_subgraph() != 0;

    // --- Exchange graph: handles all dependency-aware mesh exchanges ---
    // MeshExch4 (corrector density) is always active.
    // PreSolveExch/PostSolveExch (pressure) edges only wired when parallel pressure is on;
    // unused types sit idle until TerminationData.
    auto exchGraph = std::make_shared<ExchangeGraph<
        ExchKind<MeshState::MeshExch1, MeshState::PostPredExch>,
        ExchKind<MeshState::MeshExch2, MeshState::PostRadExch>,
        ExchKind<MeshState::MeshExch3, MeshState::PostPredVelExch>,
        ExchKind<MeshState::MeshExch4, MeshState::PostCorrStep1>,
        ExchKind<MeshState::MeshExch7, MeshState::PostParticleOps>,
        ExchKind<MeshState::PreSolveExch, MeshState::SolvePhase>,
        ExchKind<MeshState::PostSolveExch, MeshState::VelErrorPhase>>>(
        depGraph, commService, "Exchange");

    // --- Wire the graph ---

    // Input: MeshData<Init> -> TimestepTask (first iteration, emit to Predictor)
    graph->input<MeshData<MeshState::Init>>(timestepTask);
    // Input: TerminationData -> Predictor + Corrector + Exchange (for cycle termination)
    graph->input<TerminationData>(predictorSubgraph);
    graph->input<TerminationData>(correctorSubgraph);
    graph->input<TerminationData>(exchGraph);

    // Predictor -> Corrector (MeshData<> only; PredPressure routes elsewhere)
    graph->edges(predictorSubgraph, correctorSubgraph);

    // Predictor ↔ Exchange (MeshExch1 out, PostPredExch back)
    graph->edges(predictorSubgraph, exchGraph);
    graph->edges(exchGraph, predictorSubgraph);

    // Corrector ↔ Exchange (MeshExch4 out, PostCorrStep1 back)
    graph->edges(correctorSubgraph, exchGraph);
    graph->edges(exchGraph, correctorSubgraph);

    // Corrector -> TimestepTask (MeshData<>)
    graph->edges(correctorSubgraph, timestepTask);

    // Cycle: TimestepTask -> back to Predictor (MeshData<>)
    graph->edges(timestepTask, predictorSubgraph);

    // Graph output: TimestepTask emits BarrierData when simulation is done
    graph->outputs(timestepTask);

    // --- Pressure subgraph (when enabled) ---
    if (useParallelPressure) {
        auto pressureSubgraph = buildPressureIterationSubgraph(
            nmeshes, budget, fds_get_pres_flag());

        graph->input<TerminationData>(pressureSubgraph);

        // Predictor <-> PressureSubgraph (PredPressure)
        graph->edges(predictorSubgraph, pressureSubgraph);
        graph->edges(pressureSubgraph, predictorSubgraph);

        // Corrector <-> PressureSubgraph (CorrPressure)
        graph->edges(correctorSubgraph, pressureSubgraph);
        graph->edges(pressureSubgraph, correctorSubgraph);

        // PressureSubgraph <-> ExchangeGraph (PreSolveExch/PostSolveExch <-> SolvePhase/VelErrorPhase)
        graph->edges(pressureSubgraph, exchGraph);
        graph->edges(exchGraph, pressureSubgraph);
    }

    return graph;
}

#endif // FDS_GRAPH_H
