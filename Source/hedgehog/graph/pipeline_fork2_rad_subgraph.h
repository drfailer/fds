#ifndef PIPELINE_FORK2_RAD_SUBGRAPH_H
#define PIPELINE_FORK2_RAD_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/pipeline_fork2_data.h"
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "corr_radiation_subgraph.h"

/// Unwrap task: Fork2RadWork -> MeshData (for feeding into CorrRadiation sub-graph).
class Fork2RadUnwrapTask
    : public hh::AbstractTask<1, Fork2RadWork, MeshData> {
public:
    Fork2RadUnwrapTask()
        : hh::AbstractTask<1, Fork2RadWork, MeshData>("Fork2RadUnwrap", 1) {}

    void execute(std::shared_ptr<Fork2RadWork> work) override {
        this->addResult(work->meshData);
    }

    std::shared_ptr<hh::AbstractTask<1, Fork2RadWork, MeshData>>
    copy() override { return std::make_shared<Fork2RadUnwrapTask>(); }
};

/// Wrap task: BarrierData -> Fork2RadBarrier (after CorrRadiation sub-graph completes).
class Fork2RadWrapTask
    : public hh::AbstractTask<1, BarrierData, Fork2RadBarrier> {
public:
    Fork2RadWrapTask()
        : hh::AbstractTask<1, BarrierData, Fork2RadBarrier>("Fork2RadWrap", 1) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        this->addResult(std::make_shared<Fork2RadBarrier>(data));
    }

    std::shared_ptr<hh::AbstractTask<1, BarrierData, Fork2RadBarrier>>
    copy() override { return std::make_shared<Fork2RadWrapTask>(); }
};

/// Build the Fork 2 Branch C sub-graph: Fork2RadWork -> CorrRadiation -> Fork2RadBarrier.
/// Wraps the existing CorrRadiation sub-graph in type adapters for fork/join routing.
inline auto buildFork2RadSubgraph(int nmeshes, size_t kernelThreads) {
    auto subgraph = std::make_shared<hh::Graph<1, Fork2RadWork, Fork2RadBarrier>>(
        "Fork2-BranchC-Radiation");

    auto unwrap = std::make_shared<Fork2RadUnwrapTask>();
    auto corrRadSubgraph = buildCorrRadiationSubgraph(nmeshes, kernelThreads);
    auto wrap = std::make_shared<Fork2RadWrapTask>();

    subgraph->inputs(unwrap);
    subgraph->edges(unwrap, corrRadSubgraph);
    subgraph->edges(corrRadSubgraph, wrap);
    subgraph->outputs(wrap);

    return subgraph;
}

#endif // PIPELINE_FORK2_RAD_SUBGRAPH_H
