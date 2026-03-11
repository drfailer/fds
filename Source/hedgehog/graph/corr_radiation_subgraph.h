#ifndef CORR_RADIATION_SUBGRAPH_H
#define CORR_RADIATION_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../data/corr_radiation_data.h"
#include "../state/corr_radiation_state.h"
#include "../task/corr_radiation_kernel_task.h"

inline auto buildCorrRadiationSubgraph(int nmeshes, size_t kernelThreads) {
    auto subgraph = std::make_shared<
        hh::Graph<1, MeshData, MeshData>>("CorrRadiation");

    auto orchSM = std::make_shared<
        hh::StateManager<1, MeshData, CorrRadiationWork>>(
        std::make_shared<CorrRadiationOrchestrator>(nmeshes),
        "CorrRadOrch");
    auto kernelTask = std::make_shared<CorrRadiationKernelTask>(
        kernelThreads);
    auto collectorSM = std::make_shared<
        hh::StateManager<1, CorrRadiationWork, MeshData>>(
        std::make_shared<CorrRadiationCollector>(nmeshes),
        "CorrRadCollector");

    subgraph->inputs(orchSM);
    subgraph->edges(orchSM, kernelTask);
    subgraph->edges(kernelTask, collectorSM);
    subgraph->outputs(collectorSM);

    return subgraph;
}

#endif // CORR_RADIATION_SUBGRAPH_H
