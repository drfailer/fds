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
///   PredFinalOrchestrator(SYNTHETIC_TURBULENCE) -> VelocityBCEdgesTask(mesh-level) ->
///   Collector(CC_VELOCITY_BC if CC_IBM)
inline auto buildPredFinalSubgraph(int nmeshes, size_t kernelThreads) {
    using SubGraphType = hh::Graph<1, MeshData, BarrierData>;
    auto subgraph = std::make_shared<SubGraphType>("PredFinal");

    auto orchSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<PredFinalOrchestrator>(nmeshes), "PredFinalOrch");
    auto kernelTask = std::make_shared<VelocityBCEdgesTask>(kernelThreads, /*applyToEstimated=*/1);

    subgraph->inputs(orchSM);
    subgraph->edges(orchSM, kernelTask);

    if (fds_is_cc_ibm()) {
        auto collectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
            std::make_shared<PredFinalCCCollector>(nmeshes), "PredFinalCCCollector");
        subgraph->edges(kernelTask, collectorSM);
        subgraph->outputs(collectorSM);
    } else {
        auto collectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
            std::make_shared<CollectorState>(nmeshes), "PredFinalCollector");
        subgraph->edges(kernelTask, collectorSM);
        subgraph->outputs(collectorSM);
    }

    return subgraph;
}

/// Build the CorrFinal sub-graph.
///
///   VelocityBCEdgesTask(mesh-level) -> CorrFinalCollector(CC_VELOCITY_BC + outputs)
inline auto buildCorrFinalSubgraph(int nmeshes, size_t kernelThreads) {
    using SubGraphType = hh::Graph<1, MeshData, BarrierData>;
    auto subgraph = std::make_shared<SubGraphType>("CorrFinal");

    auto kernelTask = std::make_shared<VelocityBCEdgesTask>(kernelThreads, /*applyToEstimated=*/0);
    auto collectorSM = std::make_shared<hh::StateManager<1, MeshData, BarrierData>>(
        std::make_shared<CorrFinalCollector>(nmeshes), "CorrFinalCollector");

    subgraph->inputs(kernelTask);
    subgraph->edges(kernelTask, collectorSM);
    subgraph->outputs(collectorSM);

    return subgraph;
}

#endif // VELOCITY_BC_SUBGRAPH_H
