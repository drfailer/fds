#ifndef VELOCITY_BC_SUBGRAPH_H
#define VELOCITY_BC_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../task/velocity_bc_edges_task.h"
#include "../task/barrier_tasks.h"
#include "../state/barrier_state.h"
#include "../state/collector_state.h"

/// Build the PredFinal sub-graph.
///
/// PhaseTransition converted from state (PredFinalCollector) to task:
///   CollectorState (N MeshData → 1 BarrierData) + PhaseTransitionTask (BarrierData → N MeshData)
///
///   VelocityBCEdgesTask → CollectorState → PhaseTransitionTask → output
inline auto buildPredFinalSubgraph(int nmeshes, size_t kernelThreads) {
    using SubGraphType = hh::Graph<1, MeshData, MeshData>;
    auto subgraph = std::make_shared<SubGraphType>("PredFinal");

    auto kernelTask = std::make_shared<VelocityBCEdgesTask>(kernelThreads, /*applyToEstimated=*/1);
    auto collectorTask = std::make_shared<CollectorTask>(nmeshes, "PredFinalCollector");
    auto phaseTransTask = std::make_shared<PhaseTransitionTask>();

    subgraph->inputs(kernelTask);
    subgraph->edges(kernelTask, collectorTask);
    subgraph->edges(collectorTask, phaseTransTask);
    subgraph->outputs(phaseTransTask);

    return subgraph;
}

/// Build the CorrFinal sub-graph.
///
/// States converted to tasks:
///   - CorrFinalOrchestrator → CollectorState + BarrierTask (MeshExch6+CC_END_STEP)
///   - CorrFinalCollector → CollectorState + CorrFinalDumpTask (reduce+dump schedule)
///
///   CollectorState → OrchTask → VelocityBCEdgesTask → CollectorState → CorrFinalDumpTask
inline auto buildCorrFinalSubgraph(int nmeshes, size_t kernelThreads) {
    using SubGraphType = hh::Graph<1, MeshData, MeshData, BarrierData>;
    auto subgraph = std::make_shared<SubGraphType>("CorrFinal");

    bool ccIBM = fds_is_cc_ibm() != 0;

    // CorrFinalOrchestrator → CollectorTask + BarrierTask
    auto orchCollectorTask = std::make_shared<CollectorTask>(nmeshes, "CorrFinalOrchCollector");
    auto orchTask = makeBarrierTask("CorrFinalOrch",
        "CC_END_STEP\\nMESH_EXCHANGE(6)",
        [ccIBM](auto& meshes) {
            if (ccIBM) { fds_cc_end_step(meshes[0]->t, meshes[0]->dt, 0); }
            fds_mesh_exchange(6);
        });

    auto kernelTask = std::make_shared<VelocityBCEdgesTask>(
        kernelThreads, /*applyToEstimated=*/0, /*doIBEdges=*/1, /*isCorrFinal=*/true);

    // CorrFinalCollector → CollectorTask + CorrFinalDumpTask
    auto dumpCollectorTask = std::make_shared<CollectorTask>(nmeshes, "CorrFinalDumpCollector");
    auto dumpTask = std::make_shared<CorrFinalDumpTask>();

    subgraph->inputs(orchCollectorTask);
    subgraph->edges(orchCollectorTask, orchTask);
    subgraph->edges(orchTask, kernelTask);
    subgraph->edges(kernelTask, dumpCollectorTask);
    subgraph->edges(dumpCollectorTask, dumpTask);
    subgraph->outputs(dumpTask);

    return subgraph;
}

#endif // VELOCITY_BC_SUBGRAPH_H
