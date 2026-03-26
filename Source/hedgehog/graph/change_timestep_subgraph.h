#ifndef CHANGE_TIMESTEP_SUBGRAPH_H
#define CHANGE_TIMESTEP_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/change_timestep_data.h"
#include "../data/mesh_data.h"
#include "../data/termination_data.h"
#include "../task/change_timestep_tasks.h"
#include "../state/change_timestep_state.h"

/// Build the time step retry sub-graph.
///
/// Architecture with 2 tasks + loop state:
///   RetryPreKernel → RetryMomDivKernel → RetryMomDivCollector → RetryLoopSM
///                  ↘ (bypass when done=true) ──────────────────↗
///
/// RetryLoopSM includes the post-kernel work (divergence exchange, div_p2,
/// pressure iteration, velocity predictor) and the retry check.
///
/// RetryLoopSM uses type-based routing:
///   - RetrySequenceData → cycles back to RetryPreKernel for another retry
///   - MeshData → exits the sub-graph (retry complete or no retry needed)
inline auto buildChangeTimeStepSubgraph(int nmeshes, size_t kernelThreads) {
    using SubGraphType = hh::Graph<2, BarrierData, TerminationData, MeshData>;
    auto subgraph = std::make_shared<SubGraphType>("ChangeTimeStepSubgraph");

    auto retryPreKernel = std::make_shared<RetryPreKernelTask>();
    auto retryMomDivKernel = std::make_shared<RetryMomentumDivKernelTask>(
        static_cast<size_t>(nmeshes));
    auto retryMomDivCollSM = std::make_shared<hh::StateManager<
        1, MeshData, RetrySequenceData>>(
        std::make_shared<RetryMomentumDivCollector>(nmeshes), "RetryMomDivCollector");

    auto retryLoopSM = std::make_shared<RetryLoopStateManager>(
        std::make_shared<RetryLoopState>(), "RetryLoop");

    // Entry point
    subgraph->input<BarrierData>(retryPreKernel);
    subgraph->input<TerminationData>(retryLoopSM);

    // RetryPreKernel outputs:
    //   MeshData → parallel kernel (retry path)
    //   RetrySequenceData → bypass to loop state (no-retry path)
    subgraph->edges(retryPreKernel, retryMomDivKernel);      // MeshData (scatter)
    subgraph->edges(retryPreKernel, retryLoopSM);             // RetrySequenceData (bypass)

    // Kernel → collector → loop state
    subgraph->edges(retryMomDivKernel, retryMomDivCollSM);    // MeshData
    subgraph->edges(retryMomDivCollSM, retryLoopSM);           // RetrySequenceData

    // Cycle: RetryLoopState emits RetrySequenceData → back to RetryPreKernel
    subgraph->edges(retryLoopSM, retryPreKernel);

    // Exit: RetryLoopState emits MeshData → subgraph output
    subgraph->outputs(retryLoopSM);

    return subgraph;
}

#endif // CHANGE_TIMESTEP_SUBGRAPH_H
