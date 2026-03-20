#ifndef PIPELINE_FORK1_VFLUX_SUBGRAPH_H
#define PIPELINE_FORK1_VFLUX_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/pipeline_fork1_data.h"
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"
#include "../state/div_setup_state.h"
#include "../task/div_setup_kernel_task.h"
#include "velocity_flux_block_subgraph.h"

/// Unwrap task: Fork1VFluxWork -> MeshData (for feeding into VFLUX sub-graph).
class Fork1VFluxUnwrapTask
    : public hh::AbstractTask<1, Fork1VFluxWork, MeshData> {
public:
    Fork1VFluxUnwrapTask()
        : hh::AbstractTask<1, Fork1VFluxWork, MeshData>("Fork1VFluxUnwrap", 1) {}

    void execute(std::shared_ptr<Fork1VFluxWork> work) override {
        this->addResult(work->meshData);
    }

    std::shared_ptr<hh::AbstractTask<1, Fork1VFluxWork, MeshData>>
    copy() override { return std::make_shared<Fork1VFluxUnwrapTask>(); }
};

/// Wrap task: MeshData -> Fork1VFluxResult (after VFLUX sub-graph completes).
class Fork1VFluxWrapTask
    : public hh::AbstractTask<1, MeshData, Fork1VFluxResult> {
public:
    Fork1VFluxWrapTask()
        : hh::AbstractTask<1, MeshData, Fork1VFluxResult>("Fork1VFluxWrap", 1) {}

    void execute(std::shared_ptr<MeshData> data) override {
        this->addResult(std::make_shared<Fork1VFluxResult>(data));
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, Fork1VFluxResult>>
    copy() override { return std::make_shared<Fork1VFluxWrapTask>(); }
};

/// Build the Fork 1 Branch A sub-graph: Fork1VFluxWork -> VFLUX -> Fork1VFluxResult.
/// Wraps the existing VFLUX block sub-graph (or mesh-level fallback) in
/// type adapters to maintain distinct fork/join edge types.
inline auto buildFork1VFluxSubgraph(int nmeshes, size_t kernelThreads,
                                     size_t blockThreads, int numBlocks,
                                     bool canBlockFlux, bool ccIBM) {
    auto subgraph = std::make_shared<hh::Graph<1, Fork1VFluxWork, Fork1VFluxResult>>(
        "Fork1-BranchA-VFlux");

    auto unwrap = std::make_shared<Fork1VFluxUnwrapTask>();
    auto wrap = std::make_shared<Fork1VFluxWrapTask>();

    subgraph->inputs(unwrap);

    if (canBlockFlux) {
        auto vfluxBlock = buildVelocityFluxBlockSubgraph(nmeshes, blockThreads, numBlocks);
        subgraph->edges(unwrap, vfluxBlock);
        subgraph->edges(vfluxBlock, wrap);
    } else if (ccIBM) {
        // CC_IBM mesh-level: orchestrator + kernel
        auto orchSM = std::make_shared<hh::StateManager<1, MeshData, MeshData>>(
            std::make_shared<CorrDivSetupOrchestrator>(nmeshes), "Fork1VFluxOrch");
        auto kernel = std::make_shared<DivSetupKernelTask>(kernelThreads);
        subgraph->edges(unwrap, orchSM);
        subgraph->edges(orchSM, kernel);
        subgraph->edges(kernel, wrap);
    } else {
        // Mesh-level fallback (Coriolis, patch velocity, etc.)
        auto kernel = std::make_shared<DivSetupKernelTask>(kernelThreads);
        subgraph->edges(unwrap, kernel);
        subgraph->edges(kernel, wrap);
    }

    subgraph->outputs(wrap);
    return subgraph;
}

#endif // PIPELINE_FORK1_VFLUX_SUBGRAPH_H
