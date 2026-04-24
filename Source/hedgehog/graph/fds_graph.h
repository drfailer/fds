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
#include "compute_subgraph.h"
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

    auto icyc = std::make_shared<int>(1);
    auto timestepTask = std::make_shared<TimestepTask>(
        nmeshes, static_cast<size_t>(nmeshes), tEnd, icyc);

    bool useParallelPressure = fds_use_pressure_subgraph() != 0;

    auto exchGraph = std::make_shared<ExchangeGraph<
        ExchKind<MeshState::MeshExch1, MeshState::PostPredExch>,
        ExchKind<MeshState::MeshExch2, MeshState::PostRadExch>,
        ExchKind<MeshState::MeshExch3, MeshState::PostPredVelExch>,
        ExchKind<MeshState::MeshExch4, MeshState::PostCorrStep1>,
        ExchKind<MeshState::MeshExch7, MeshState::PostParticleOps>,
        ExchKind<MeshState::PreSolveExch, MeshState::SolvePhase>,
        ExchKind<MeshState::PostSolveExch, MeshState::VelErrorPhase>>>(
        depGraph, commService, "Exchange");

    graph->input<MeshData<MeshState::Init>>(timestepTask);
    graph->input<TerminationData>(exchGraph);
    graph->outputs(timestepTask);

    bool ccIBM = fds_is_cc_ibm() != 0;

    if (!ccIBM) {
        // --- Two-lane compute subgraph (non-CC_IBM) ---
        auto computeSubgraph = buildComputeSubgraph(nmeshes, budget);

        graph->input<TerminationData>(computeSubgraph);

        // TimestepTask ↔ Compute (MeshData<> cycle)
        graph->edges(timestepTask, computeSubgraph);
        graph->edges(computeSubgraph, timestepTask);

        // Compute ↔ Exchange (MeshExch1-7 out, Post*Exch back)
        graph->edges(computeSubgraph, exchGraph);
        graph->edges(exchGraph, computeSubgraph);

        if (useParallelPressure) {
            auto pressureSubgraph = buildPressureIterationSubgraph(
                nmeshes, budget, fds_get_pres_flag());
            graph->input<TerminationData>(pressureSubgraph);
            graph->edges(computeSubgraph, pressureSubgraph);
            graph->edges(pressureSubgraph, computeSubgraph);
            graph->edges(pressureSubgraph, exchGraph);
            graph->edges(exchGraph, pressureSubgraph);
        }
    } else {
        // --- CC_IBM: separate predictor + corrector subgraphs ---
        auto predictorSubgraph = buildPredictorSubgraph(nmeshes, budget);
        auto correctorSubgraph = buildCorrectorSubgraph(nmeshes, budget);

        graph->input<TerminationData>(predictorSubgraph);
        graph->input<TerminationData>(correctorSubgraph);

        graph->edges(predictorSubgraph, correctorSubgraph);
        graph->edges(predictorSubgraph, exchGraph);
        graph->edges(exchGraph, predictorSubgraph);
        graph->edges(correctorSubgraph, exchGraph);
        graph->edges(exchGraph, correctorSubgraph);
        graph->edges(correctorSubgraph, timestepTask);
        graph->edges(timestepTask, predictorSubgraph);

        if (useParallelPressure) {
            auto pressureSubgraph = buildPressureIterationSubgraph(
                nmeshes, budget, fds_get_pres_flag());
            graph->input<TerminationData>(pressureSubgraph);
            graph->edges(predictorSubgraph, pressureSubgraph);
            graph->edges(pressureSubgraph, predictorSubgraph);
            graph->edges(correctorSubgraph, pressureSubgraph);
            graph->edges(pressureSubgraph, correctorSubgraph);
            graph->edges(pressureSubgraph, exchGraph);
            graph->edges(exchGraph, pressureSubgraph);
        }
    }

    return graph;
}

#endif // FDS_GRAPH_H
