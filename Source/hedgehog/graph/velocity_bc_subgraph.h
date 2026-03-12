#ifndef VELOCITY_BC_SUBGRAPH_H
#define VELOCITY_BC_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../task/velocity_bc_edges_task.h"
#include "../state/velocity_bc_state.h"

/// Build the PredFinal sub-graph.
///
/// Architecture:
///   1. Sequential preprocessing (PredFinalOrchestrator):
///      - SYNTHETIC_TURBULENCE_IF_ENABLED (SEM inflow BC — uses RANDOM_NUMBER)
///
///   2. Parallel kernel execution (VelocityBCEdgesTask, applyToEstimated=1):
///      - MATCH_VELOCITY_KERNEL, VELOCITY_BC_PREPROCESSING, VELOCITY_BC_PROCESS_EDGES_KERNEL
///
///   3. CC_IBM only: Sequential finalization (PredFinalCCCollector):
///      - CC_VELOCITY_BC (cut-cell velocity BC)
inline auto buildPredFinalSubgraph(int nmeshes, size_t kernelThreads) {
    using SubGraphType = hh::Graph<1, MeshData, MeshData>;
    auto subgraph = std::make_shared<SubGraphType>("PredFinal");

    auto orchSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<PredFinalOrchestrator>(nmeshes), "PredFinalOrch");
    auto kernelTask = std::make_shared<VelocityBCEdgesTask>(kernelThreads, /*applyToEstimated=*/1);

    subgraph->inputs(orchSM);
    subgraph->edges(orchSM, kernelTask);

    if (fds_is_cc_ibm()) {
        auto collectorSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
            std::make_shared<PredFinalCCCollector>(nmeshes), "PredFinalCCCollector");
        subgraph->edges(kernelTask, collectorSM);
        subgraph->outputs(collectorSM);
    } else {
        subgraph->outputs(kernelTask);
    }

    return subgraph;
}

/// Build the CorrFinal sub-graph.
///
/// Architecture:
///   1. Parallel kernel execution (VelocityBCEdgesTask, applyToEstimated=0):
///      - MATCH_VELOCITY_KERNEL, VELOCITY_BC_PREPROCESSING, VELOCITY_BC_PROCESS_EDGES_KERNEL
///
///   2. Sequential finalization (CorrFinalCollector):
///      - CC_VELOCITY_BC (if CC_IBM active — no-op otherwise)
///      - UPDATE_GLOBAL_OUTPUTS (per-mesh output accumulation)
inline auto buildCorrFinalSubgraph(int nmeshes, size_t kernelThreads) {
    using SubGraphType = hh::Graph<1, MeshData, MeshData>;
    auto subgraph = std::make_shared<SubGraphType>("CorrFinal");

    auto kernelTask = std::make_shared<VelocityBCEdgesTask>(kernelThreads, /*applyToEstimated=*/0);
    auto collectorSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
        std::make_shared<CorrFinalCollector>(nmeshes), "CorrFinalCollector");

    subgraph->inputs(kernelTask);
    subgraph->edges(kernelTask, collectorSM);
    subgraph->outputs(collectorSM);

    return subgraph;
}

#endif // VELOCITY_BC_SUBGRAPH_H
