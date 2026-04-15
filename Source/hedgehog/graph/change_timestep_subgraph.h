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
/// RetryLoopSM collects N MeshData<> from the kernel, runs post-kernel work
/// (divergence exchange, div_p2, pressure iteration, velocity predictor),
/// checks for retry, and on exit runs CC_END_STEP + MESH_EXCHANGE(3).
///
/// Type-specific edges avoid routing MeshData<> from RetryPreKernel to
/// RetryLoopSM (only RetrySequenceData should take that path).
///
/// Output types (type-based routing):
///   - RetrySequenceData → cycles back to RetryPreKernel for another retry
///   - MeshData<> → exits the sub-graph (retry complete or no retry needed)
inline auto buildChangeTimeStepSubgraph(int nmeshes, size_t kernelThreads, bool ccIBM) {
    using SubGraphType = hh::Graph<2, MeshData<>, TerminationData, MeshData<>>;
    auto subgraph = std::make_shared<SubGraphType>("ChangeTimeStepSubgraph");

    auto retryPreKernel = std::make_shared<RetryPreKernelTask>(nmeshes);
    auto retryMomDivKernel = std::make_shared<RetryMomentumDivKernelTask>(
        kernelThreads);

    auto retryLoopSM = std::make_shared<RetryLoopStateManager>(
        std::make_shared<RetryLoopState>(nmeshes, ccIBM), "RetryLoop");

    // Entry point: MeshData<> from VelocityPredictor → RetryPreKernel (collects N)
    subgraph->input<MeshData<>>(retryPreKernel);
    subgraph->input<TerminationData>(retryLoopSM);

    // RetryPreKernel → kernel: MeshData<> scatter for parallel processing
    // Must use edge<MeshData<>> to avoid routing RetrySequenceData to kernel
    subgraph->template edge<MeshData<>>(retryPreKernel, retryMomDivKernel);

    // RetryPreKernel → loop: RetrySequenceData ONLY (bypass path)
    subgraph->template edge<RetrySequenceData>(retryPreKernel, retryLoopSM);

    // Kernel → loop: MeshData<> (collected internally by merged state)
    subgraph->edges(retryMomDivKernel, retryLoopSM);

    // Cycle: loop → pre-kernel: RetrySequenceData ONLY
    // Must use edge<T> to avoid routing MeshData<> back to pre-kernel
    subgraph->template edge<RetrySequenceData>(retryLoopSM, retryPreKernel);

    // Exit: RetryLoopState emits MeshData<> → subgraph output
    subgraph->outputs(retryLoopSM);

    return subgraph;
}

#endif // CHANGE_TIMESTEP_SUBGRAPH_H
