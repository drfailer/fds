#ifndef VELOCITY_BC_SUBGRAPH_H
#define VELOCITY_BC_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/velocity_bc_data.h"
#include "../task/velocity_bc_edges_task.h"
#include "../state/velocity_bc_state.h"

/// Build the PredFinal sub-graph (Pattern B: complex routine parallelization).
///
/// Three-phase architecture:
///   1. Sequential preprocessing (PredFinalOrchestrator):
///      - MATCH_VELOCITY (cross-mesh velocity interpolation)
///      - SYNTHETIC_TURBULENCE_IF_ENABLED (SEM inflow BC)
///      - VELOCITY_BC_PREPROCESSING (OMESH reads for wall boundary velocities)
///
///   2. Parallel kernel execution (VelocityBCEdgesTask):
///      - VELOCITY_BC_PROCESS_EDGES_KERNEL (all edge boundary conditions)
///      - Thread-safe: uses explicit M% access
///
///   3. Sequential finalization (PredFinalCollector):
///      - CC_VELOCITY_BC (cut-cell velocity BC if CC_IBM active)
///
/// Replaces the sequential PredFinalTask.
inline auto buildPredFinalSubgraph(int nmeshes, size_t kernelThreads) {
    using SubGraphType = hh::Graph<1, MeshData, MeshData>;
    auto subgraph = std::make_shared<SubGraphType>("PredFinal");

    auto orchSM = std::make_shared<hh::StateManager<1, MeshData, VelocityBCWork>>(
        std::make_shared<PredFinalOrchestrator>(nmeshes), "PredFinalOrch");
    auto kernelTask = std::make_shared<VelocityBCEdgesTask>(kernelThreads);
    auto collectorSM = std::make_shared<hh::StateManager<1, VelocityBCWork, MeshData>>(
        std::make_shared<PredFinalCollector>(nmeshes), "PredFinalCollector");

    subgraph->inputs(orchSM);
    subgraph->edges(orchSM, kernelTask);
    subgraph->edges(kernelTask, collectorSM);
    subgraph->outputs(collectorSM);

    return subgraph;
}

/// Build the CorrFinal sub-graph (Pattern B: complex routine parallelization).
///
/// Three-phase architecture:
///   1. Sequential preprocessing (CorrFinalOrchestrator):
///      - MATCH_VELOCITY (cross-mesh velocity interpolation)
///      - VELOCITY_BC_PREPROCESSING (OMESH reads for wall boundary velocities)
///
///   2. Parallel kernel execution (VelocityBCEdgesTask):
///      - VELOCITY_BC_PROCESS_EDGES_KERNEL (all edge boundary conditions)
///      - Thread-safe: uses explicit M% access
///
///   3. Sequential finalization (CorrFinalCollector):
///      - CC_VELOCITY_BC (cut-cell velocity BC if CC_IBM active)
///      - UPDATE_GLOBAL_OUTPUTS (per-mesh output accumulation)
///
/// Replaces the sequential CorrFinalTask.
inline auto buildCorrFinalSubgraph(int nmeshes, size_t kernelThreads) {
    using SubGraphType = hh::Graph<1, MeshData, MeshData>;
    auto subgraph = std::make_shared<SubGraphType>("CorrFinal");

    auto orchSM = std::make_shared<hh::StateManager<1, MeshData, VelocityBCWork>>(
        std::make_shared<CorrFinalOrchestrator>(nmeshes), "CorrFinalOrch");
    auto kernelTask = std::make_shared<VelocityBCEdgesTask>(kernelThreads);
    auto collectorSM = std::make_shared<hh::StateManager<1, VelocityBCWork, MeshData>>(
        std::make_shared<CorrFinalCollector>(nmeshes), "CorrFinalCollector");

    subgraph->inputs(orchSM);
    subgraph->edges(orchSM, kernelTask);
    subgraph->edges(kernelTask, collectorSM);
    subgraph->outputs(collectorSM);

    return subgraph;
}

#endif // VELOCITY_BC_SUBGRAPH_H
