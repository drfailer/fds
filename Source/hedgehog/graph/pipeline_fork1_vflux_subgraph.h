#ifndef PIPELINE_FORK1_VFLUX_SUBGRAPH_H
#define PIPELINE_FORK1_VFLUX_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"
#include "../state/div_setup_state.h"
#include "../task/div_setup_kernel_task.h"

/// Build the Fork 1 Branch A sub-graph: MeshData<> -> VFLUX -> MeshData<>.
inline auto buildFork1VFluxSubgraph(int nmeshes, bool ccIBM,
                                     size_t divSetupThreads) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData<>, MeshData<>>>(
        "Fork1-BranchA-VFlux");

    if (ccIBM) {
        auto orchTask = std::make_shared<CorrDivSetupOrchestrator>(nmeshes);
        auto kernel = std::make_shared<DivSetupKernelTask>(divSetupThreads);
        subgraph->inputs(orchTask);
        subgraph->edges(orchTask, kernel);
        subgraph->outputs(kernel);
    } else {
        auto kernel = std::make_shared<DivSetupKernelTask>(divSetupThreads);
        subgraph->inputs(kernel);
        subgraph->outputs(kernel);
    }

    return subgraph;
}

#endif // PIPELINE_FORK1_VFLUX_SUBGRAPH_H
