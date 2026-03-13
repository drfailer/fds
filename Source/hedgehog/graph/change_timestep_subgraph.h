#ifndef CHANGE_TIMESTEP_SUBGRAPH_H
#define CHANGE_TIMESTEP_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/change_timestep_data.h"
#include "../data/mesh_data.h"
#include "../task/change_timestep_tasks.h"
#include "../state/change_timestep_state.h"

/// Build the time step retry sub-graph.
///
/// Consolidated architecture with 3 tasks + loop state:
///   RetryPreKernel → RetryMomDivKernel → RetryMomDivCollector → RetryPostKernel → RetryLoopSM
///                  ↘ (bypass when done=true) ──────────────────↗
///
/// RetryPreKernelTask accepts BarrierData (first entry) and RetrySequenceData (cycle).
/// When no retry needed: emits RetrySequenceData(done=true) → bypasses kernel → RetryPostKernel passes through.
/// When retry needed: runs all sequential pre-kernel work, scatters MeshData for parallel kernel.
///
/// RetryLoopSM uses type-based routing:
///   - RetrySequenceData → cycles back to RetryPreKernel for another retry
///   - MeshData → exits the sub-graph (retry complete or no retry needed)
///
/// @param tEnd Simulation end time (used for data-driven cycle termination)
/// @param nmeshes Number of local meshes
/// @param kernelThreads Number of threads for parallel kernel tasks
/// @return Shared pointer to the constructed sub-graph
inline auto buildChangeTimeStepSubgraph(double tEnd, int nmeshes,
                                         size_t kernelThreads) {
    using SubGraphType = hh::Graph<1, BarrierData, MeshData>;
    auto subgraph = std::make_shared<SubGraphType>("ChangeTimeStepSubgraph");

    // --- Create tasks ---
    auto retryPreKernel = std::make_shared<RetryPreKernelTask>();
    auto retryMomDivKernel = std::make_shared<RetryMomentumDivKernelTask>(kernelThreads);
    auto retryMomDivCollSM = std::make_shared<hh::StateManager<
        1, MeshData, RetrySequenceData>>(
        std::make_shared<RetryMomentumDivCollector>(nmeshes), "RetryMomDivCollector");
    auto retryPostKernel = std::make_shared<RetryPostKernelTask>();

    // --- Create retry loop state manager (data-driven canTerminate) ---
    auto retryLoopSM = std::make_shared<RetryLoopStateManager>(
        std::make_shared<RetryLoopState>(tEnd), "RetryLoop");

    // --- Wire the sub-graph ---

    // Entry point: RetryPreKernelTask receives BarrierData
    subgraph->inputs(retryPreKernel);

    // RetryPreKernel outputs:
    //   MeshData → parallel kernel (retry path)
    //   RetrySequenceData → bypass to post-kernel (no-retry path)
    subgraph->edges(retryPreKernel, retryMomDivKernel);      // MeshData (scatter)
    subgraph->edges(retryPreKernel, retryPostKernel);         // RetrySequenceData (bypass)

    // Kernel → collector → post-kernel
    subgraph->edges(retryMomDivKernel, retryMomDivCollSM);    // MeshData
    subgraph->edges(retryMomDivCollSM, retryPostKernel);      // RetrySequenceData

    // Post-kernel → retry loop state
    subgraph->edges(retryPostKernel, retryLoopSM);            // RetrySequenceData

    // Cycle: RetryLoopState emits RetrySequenceData → back to RetryPreKernel
    subgraph->edges(retryLoopSM, retryPreKernel);

    // Exit: RetryLoopState emits MeshData → subgraph output
    subgraph->outputs(retryLoopSM);

    return subgraph;
}

#endif // CHANGE_TIMESTEP_SUBGRAPH_H
