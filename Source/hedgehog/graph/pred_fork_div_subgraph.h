#ifndef PRED_FORK_DIV_SUBGRAPH_H
#define PRED_FORK_DIV_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../task/pred_fork_tasks.h"
#include "wallbc_subgraph.h"
#include "wallbc_block_subgraph.h"

/// Build the Predictor Fork Branch B sub-graph: WALL_BC -> DIV_P1_early.
/// Only used for non-CC_IBM.
inline auto buildPredForkDivSubgraph(int nmeshes, size_t blockThreads,
                                      int numBlocks, bool canBlockWallBC) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData, MeshData>>(
        "PredFork-BranchB-WallBC+DivEarly");

    auto divEarly = std::make_shared<DivP1EarlyTask>(static_cast<size_t>(nmeshes));

    auto wallBC = canBlockWallBC
        ? buildWallBCBlockSubgraph(nmeshes, blockThreads, numBlocks)
        : buildWallBCSubgraph(nmeshes, static_cast<size_t>(nmeshes));

    subgraph->inputs(wallBC);
    subgraph->edges(wallBC, divEarly);
    subgraph->outputs(divEarly);
    return subgraph;
}

#endif // PRED_FORK_DIV_SUBGRAPH_H
