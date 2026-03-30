#ifndef VELOCITY_BC_SUBGRAPH_H
#define VELOCITY_BC_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../task/velocity_bc_edges_task.h"
#include "../state/velocity_bc_state.h"
#include "../state/collector_state.h"

/// Build the PredFinal sub-graph.
///
/// MeshExchange(3) + SyntheticTurbulence merged into upstream barrier.
/// PhaseTransition merged into the collector.
///
///   VelocityBCEdgesTask(mesh-level) -> PredFinalCollector(phase transition)
///
/// Outputs MeshData directly (no BarrierData intermediate).
inline auto buildPredFinalSubgraph(int nmeshes, size_t kernelThreads) {
    using SubGraphType = hh::Graph<1, MeshData, MeshData>;
    auto subgraph = std::make_shared<SubGraphType>("PredFinal");

    auto kernelTask = std::make_shared<VelocityBCEdgesTask>(kernelThreads, /*applyToEstimated=*/1);
    auto collectorSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<PredFinalCollector>(nmeshes), "PredFinalCollector");

    subgraph->inputs(kernelTask);
    subgraph->edges(kernelTask, collectorSM);
    subgraph->outputs(collectorSM);

    return subgraph;
}

/// Build the CorrFinal sub-graph.
///
/// MeshExchange(6b) merged into the orchestrator.
///
///   CorrFinalOrch(MeshExch6+CC_END_STEP) -> VelocityBCEdgesTask -> CorrFinalCollector
inline auto buildCorrFinalSubgraph(int nmeshes, size_t kernelThreads) {
    using SubGraphType = hh::Graph<1, MeshData, BarrierData>;
    auto subgraph = std::make_shared<SubGraphType>("CorrFinal");

    bool ccIBM = fds_is_cc_ibm() != 0;
    auto orchSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<CorrFinalOrchestrator>(nmeshes, ccIBM), "CorrFinalOrch");
    auto kernelTask = std::make_shared<VelocityBCEdgesTask>(kernelThreads, /*applyToEstimated=*/0);
    auto collectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CorrFinalCollector>(nmeshes), "CorrFinalCollector");

    subgraph->inputs(orchSM);
    subgraph->edges(orchSM, kernelTask);
    subgraph->edges(kernelTask, collectorSM);
    subgraph->outputs(collectorSM);

    return subgraph;
}

#endif // VELOCITY_BC_SUBGRAPH_H
