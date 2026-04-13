#ifndef CORR_RADIATION_SUBGRAPH_H
#define CORR_RADIATION_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../data/corr_radiation_data.h"
#include "../state/corr_radiation_state.h"
#include "../task/corr_radiation_kernel_task.h"

inline auto buildCorrRadiationSubgraph(int nmeshes, size_t kernelThreads) {
    auto subgraph = std::make_shared<
        hh::Graph<1, MeshData<>, MeshData<>>>("CorrRadiation");

    auto orchTask = std::make_shared<CorrRadiationOrchestrator>(nmeshes);
    auto kernelTask = std::make_shared<CorrRadiationKernelTask>(
        kernelThreads);
    auto collectorTask = std::make_shared<CorrRadiationCollector>(nmeshes);

    subgraph->inputs(orchTask);
    subgraph->edges(orchTask, kernelTask);
    subgraph->edges(kernelTask, collectorTask);
    subgraph->outputs(collectorTask);

    return subgraph;
}

#endif // CORR_RADIATION_SUBGRAPH_H
