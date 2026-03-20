#ifndef PRED_FORK_DIV_SUBGRAPH_H
#define PRED_FORK_DIV_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/pred_fork_data.h"
#include "../data/mesh_data.h"
#include "../task/pred_fork_tasks.h"
#include "wallbc_subgraph.h"
#include "wallbc_block_subgraph.h"

/// Unwrap task: PredForkDivWork -> MeshData.
class PredForkDivUnwrapTask
    : public hh::AbstractTask<1, PredForkDivWork, MeshData> {
public:
    PredForkDivUnwrapTask()
        : hh::AbstractTask<1, PredForkDivWork, MeshData>(
              "PredForkDivUnwrap", 1) {}

    void execute(std::shared_ptr<PredForkDivWork> work) override {
        this->addResult(work->meshData);
    }

    std::shared_ptr<hh::AbstractTask<1, PredForkDivWork, MeshData>>
    copy() override { return std::make_shared<PredForkDivUnwrapTask>(); }
};

/// Wrap task: MeshData -> PredForkDivResult.
class PredForkDivWrapTask
    : public hh::AbstractTask<1, MeshData, PredForkDivResult> {
public:
    PredForkDivWrapTask()
        : hh::AbstractTask<1, MeshData, PredForkDivResult>(
              "PredForkDivWrap", 1) {}

    void execute(std::shared_ptr<MeshData> data) override {
        this->addResult(std::make_shared<PredForkDivResult>(data));
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, PredForkDivResult>>
    copy() override { return std::make_shared<PredForkDivWrapTask>(); }
};

/// Build the Predictor Fork Branch B sub-graph: WALL_BC -> DIV_P1_early.
/// Wraps the existing WallBC sub-graph (block or mesh-level) in type
/// adapters, then chains DIV_P1 early phases (PHASE=2, WORK_BRANCH=2).
///
/// Only used for non-CC_IBM (CC_IBM falls back to sequential predictor path).
inline auto buildPredForkDivSubgraph(int nmeshes, size_t kernelThreads,
                                      size_t blockThreads, int numBlocks,
                                      bool canBlockWallBC) {
    auto subgraph = std::make_shared<hh::Graph<1, PredForkDivWork, PredForkDivResult>>(
        "PredFork-BranchB-WallBC+DivEarly");

    auto unwrap = std::make_shared<PredForkDivUnwrapTask>();
    auto divEarly = std::make_shared<DivP1EarlyTask>(kernelThreads);
    auto wrap = std::make_shared<PredForkDivWrapTask>();

    auto wallBC = canBlockWallBC
        ? buildWallBCBlockSubgraph(nmeshes, blockThreads, numBlocks)
        : buildWallBCSubgraph(nmeshes, kernelThreads);

    subgraph->inputs(unwrap);
    subgraph->edges(unwrap, wallBC);
    subgraph->edges(wallBC, divEarly);
    subgraph->edges(divEarly, wrap);
    subgraph->outputs(wrap);
    return subgraph;
}

#endif // PRED_FORK_DIV_SUBGRAPH_H
