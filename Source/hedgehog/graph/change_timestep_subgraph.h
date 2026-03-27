#ifndef CHANGE_TIMESTEP_SUBGRAPH_H
#define CHANGE_TIMESTEP_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/barrier_data.h"
#include "../data/change_timestep_data.h"
#include "../data/mesh_data.h"
#include "../data/termination_data.h"
#include "../task/change_timestep_tasks.h"
#include "../state/change_timestep_state.h"

/// Build the time step retry sub-graph.
///
/// Architecture with 2 tasks + merged collector/loop state:
///   RetryPreKernel → RetryMomDivKernel → RetryLoopSM (collects N + post-kernel)
///                  ↘ (bypass when done=true) ↗
///
/// RetryLoopSM collects N MeshData from the kernel, runs post-kernel work
/// (divergence exchange, div_p2, pressure iteration, velocity predictor),
/// and checks for retry.
///
/// Type-specific edges avoid routing MeshData from RetryPreKernel to
/// RetryLoopSM (only RetrySequenceData should take that path).
///
/// Output types (type-based routing):
///   - RetrySequenceData → cycles back to RetryPreKernel for another retry
///   - BarrierData → exits the sub-graph (retry complete or no retry needed)
inline auto buildChangeTimeStepSubgraph(int nmeshes, size_t kernelThreads) {
    using SubGraphType = hh::Graph<2, BarrierData, TerminationData, BarrierData>;
    auto subgraph = std::make_shared<SubGraphType>("ChangeTimeStepSubgraph");

    auto retryPreKernel = std::make_shared<RetryPreKernelTask>();
    auto retryMomDivKernel = std::make_shared<RetryMomentumDivKernelTask>(
        static_cast<size_t>(nmeshes));

    auto retryLoopSM = std::make_shared<RetryLoopStateManager>(
        std::make_shared<RetryLoopState>(nmeshes), "RetryLoop");

    // Entry point
    subgraph->input<BarrierData>(retryPreKernel);
    subgraph->input<TerminationData>(retryLoopSM);

    // RetryPreKernel → kernel: MeshData (all matching types, kernel only accepts MeshData)
    subgraph->edges(retryPreKernel, retryMomDivKernel);

    // RetryPreKernel → loop: RetrySequenceData ONLY (bypass path)
    // Must use edge<T> to avoid routing MeshData to the loop state
    subgraph->template edge<RetrySequenceData>(retryPreKernel, retryLoopSM);

    // Kernel → loop: MeshData (collected internally by merged state)
    subgraph->edges(retryMomDivKernel, retryLoopSM);

    // Cycle: loop → pre-kernel: RetrySequenceData ONLY
    // Must use edge<T> to avoid routing BarrierData back to pre-kernel
    subgraph->template edge<RetrySequenceData>(retryLoopSM, retryPreKernel);

    // Exit: RetryLoopState emits BarrierData → subgraph output
    subgraph->outputs(retryLoopSM);

    return subgraph;
}

#endif // CHANGE_TIMESTEP_SUBGRAPH_H
