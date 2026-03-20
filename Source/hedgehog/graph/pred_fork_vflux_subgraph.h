#ifndef PRED_FORK_VFLUX_SUBGRAPH_H
#define PRED_FORK_VFLUX_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/pred_fork_data.h"
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"
#include "../task/div_setup_kernel_task.h"
#include "../task/pred_fork_tasks.h"
#include "velocity_flux_block_subgraph.h"

/// Unwrap task: PredForkVFluxWork -> MeshData.
class PredForkVFluxUnwrapTask
    : public hh::AbstractTask<1, PredForkVFluxWork, MeshData> {
public:
    PredForkVFluxUnwrapTask()
        : hh::AbstractTask<1, PredForkVFluxWork, MeshData>(
              "PredForkVFluxUnwrap", 1) {}

    void execute(std::shared_ptr<PredForkVFluxWork> work) override {
        this->addResult(work->meshData);
    }

    std::shared_ptr<hh::AbstractTask<1, PredForkVFluxWork, MeshData>>
    copy() override { return std::make_shared<PredForkVFluxUnwrapTask>(); }
};

/// Wrap task: MeshData -> PredForkVFluxResult.
class PredForkVFluxWrapTask
    : public hh::AbstractTask<1, MeshData, PredForkVFluxResult> {
public:
    PredForkVFluxWrapTask()
        : hh::AbstractTask<1, MeshData, PredForkVFluxResult>(
              "PredForkVFluxWrap", 1) {}

    void execute(std::shared_ptr<MeshData> data) override {
        this->addResult(std::make_shared<PredForkVFluxResult>(data));
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, PredForkVFluxResult>>
    copy() override { return std::make_shared<PredForkVFluxWrapTask>(); }
};

/// Build the Predictor Fork Branch A sub-graph: VFLUX -> PARTICLE_MOMENTUM.
/// Wraps the existing VFLUX block sub-graph (or mesh-level fallback) in
/// type adapters, then chains PARTICLE_MOMENTUM at mesh level.
///
/// Only used for non-CC_IBM (CC_IBM falls back to sequential predictor path).
inline auto buildPredForkVFluxSubgraph(int nmeshes, size_t kernelThreads,
                                        size_t blockThreads, int numBlocks,
                                        bool canBlockFlux) {
    auto subgraph = std::make_shared<hh::Graph<1, PredForkVFluxWork, PredForkVFluxResult>>(
        "PredFork-BranchA-VFlux+PMom");

    auto unwrap = std::make_shared<PredForkVFluxUnwrapTask>();
    auto partMom = std::make_shared<PredPartMomKernelTask>(kernelThreads);
    auto wrap = std::make_shared<PredForkVFluxWrapTask>();

    subgraph->inputs(unwrap);

    if (canBlockFlux) {
        auto vfluxBlock = buildVelocityFluxBlockSubgraph(
            nmeshes, blockThreads, numBlocks);
        subgraph->edges(unwrap, vfluxBlock);
        subgraph->edges(vfluxBlock, partMom);
    } else {
        auto kernel = std::make_shared<DivSetupKernelTask>(kernelThreads);
        subgraph->edges(unwrap, kernel);
        subgraph->edges(kernel, partMom);
    }

    subgraph->edges(partMom, wrap);
    subgraph->outputs(wrap);
    return subgraph;
}

#endif // PRED_FORK_VFLUX_SUBGRAPH_H
