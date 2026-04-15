#ifndef PIPELINE_FORK1_VFLUX_SUBGRAPH_H
#define PIPELINE_FORK1_VFLUX_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../task/div_setup_kernel_task.h"

/// Build the Fork 1 Branch A sub-graph: MeshData<> -> VFLUX -> MeshData<>.
inline auto buildFork1VFluxSubgraph(int nmeshes, bool /*ccIBM*/,
                                     size_t divSetupThreads) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData<>, MeshData<>>>(
        "Fork1-BranchA-VFlux");

    // CC_IBM and non-CC_IBM both go directly to kernel (no orchestrator needed —
    // CC_VELOCITY_BC already moved to thread-safe DivSetupKernelTask)
    auto kernel = std::make_shared<DivSetupKernelTask>(divSetupThreads);
    subgraph->inputs(kernel);
    subgraph->outputs(kernel);

    return subgraph;
}

#endif // PIPELINE_FORK1_VFLUX_SUBGRAPH_H
