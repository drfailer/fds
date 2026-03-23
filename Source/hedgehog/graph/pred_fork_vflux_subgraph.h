#ifndef PRED_FORK_VFLUX_SUBGRAPH_H
#define PRED_FORK_VFLUX_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"
#include "../task/div_setup_kernel_task.h"
#include "../task/pred_fork_tasks.h"
#include "velocity_flux_block_subgraph.h"

/// Build the Predictor Fork Branch A sub-graph: VFLUX -> PARTICLE_MOMENTUM.
/// Only used for non-CC_IBM.
inline auto buildPredForkVFluxSubgraph(int nmeshes, size_t blockThreads,
                                        int numBlocks, bool canBlockFlux) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData, MeshData>>(
        "PredFork-BranchA-VFlux+PMom");

    auto partMom = std::make_shared<PredPartMomKernelTask>(static_cast<size_t>(nmeshes));

    if (canBlockFlux) {
        auto vfluxBlock = buildVelocityFluxBlockSubgraph(
            nmeshes, blockThreads, numBlocks);
        subgraph->inputs(vfluxBlock);
        subgraph->edges(vfluxBlock, partMom);
    } else {
        auto kernel = std::make_shared<DivSetupKernelTask>(static_cast<size_t>(nmeshes));
        subgraph->inputs(kernel);
        subgraph->edges(kernel, partMom);
    }

    subgraph->outputs(partMom);
    return subgraph;
}

#endif // PRED_FORK_VFLUX_SUBGRAPH_H
