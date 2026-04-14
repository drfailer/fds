#ifndef PRED_FORK_DIV_SUBGRAPH_H
#define PRED_FORK_DIV_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../task/pred_fork_tasks.h"
#include "../task/wallbc_kernel_task.h"

/// Build the Predictor Fork Branch B sub-graph: WALL_BC -> DIV_P1_early.
/// Only used for non-CC_IBM.
///
/// WallBC inlined: no orchestrator barrier needed in predictor (dt_bc=0, call_ht_1d=0
/// are MeshData<> defaults). WallBCKernelTask includes finalize (all per-mesh, thread-safe).
inline auto buildPredForkDivSubgraph(int nmeshes,
                                      size_t wallBCThreads,
                                      size_t divP1EarlyThreads) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData<>, MeshData<>>>(
        "PredFork-BranchB-WallBC+DivEarly");

    auto wallBCKernel = std::make_shared<WallBCKernelTask>(wallBCThreads);
    auto divEarly = std::make_shared<DivP1EarlyTask>(divP1EarlyThreads);

    subgraph->inputs(wallBCKernel);
    subgraph->edges(wallBCKernel, divEarly);
    subgraph->outputs(divEarly);
    return subgraph;
}

#endif // PRED_FORK_DIV_SUBGRAPH_H
