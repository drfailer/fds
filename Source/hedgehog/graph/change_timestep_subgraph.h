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
/// This sub-graph replaces the monolithic ChangeTimeStepTask with a clean
/// dataflow architecture:
///   1. CheckRetryTask: determines if retry is needed
///   2. Retry sequence: small focused tasks for each operation
///   3. RetryLoopState: manages the loop (cycle or exit)
///
/// The retry loop re-runs the predictor sequence with reduced time step
/// until CFL compliance is achieved.
///
/// @return Shared pointer to the constructed sub-graph
inline auto buildChangeTimeStepSubgraph() {
    using SubGraphType = hh::Graph<1, BarrierData, MeshData>;
    auto subgraph = std::make_shared<SubGraphType>("ChangeTimeStepSubgraph");

    // --- Create tasks for retry sequence ---
    auto checkRetry = std::make_shared<CheckRetryTask>();
    auto retryDensity = std::make_shared<RetryDensityTask>();
    auto retryCCDensity = std::make_shared<RetryCCDensityTask>();
    auto retryVelocityFlux = std::make_shared<RetryVelocityFluxTask>();
    auto retryHvac = std::make_shared<RetryHvacTask>();
    auto retryInitDiv = std::make_shared<RetryInitDivTask>();
    auto retryDivPart1 = std::make_shared<RetryDivergencePart1Task>();
    auto retryDivExchange = std::make_shared<RetryDivExchangeTask>();
    auto retryDivPart2 = std::make_shared<RetryDivergencePart2Task>();
    auto retryPressure = std::make_shared<RetryPressureTask>();
    auto retryVelocityPredictor = std::make_shared<RetryVelocityPredictorTask>();
    auto retryExit = std::make_shared<RetryExitTask>();

    // --- Create retry loop state manager ---
    auto retryLoopSM = std::make_shared<hh::StateManager<1, RetrySequenceData, RetrySequenceData>>(
        std::make_shared<RetryLoopState>(), "RetryLoop");

    // --- Wire the sub-graph ---

    // Entry point: CheckRetryTask receives BarrierData
    subgraph->inputs(checkRetry);

    // CheckRetryTask outputs RetrySequenceData with two paths:
    // 1. done=false → retry sequence (retryDensity will process)
    // 2. done=true → direct to exit (retryExit will process)
    subgraph->edges(checkRetry, retryDensity);  // Main path
    subgraph->edges(checkRetry, retryExit);     // Bypass path for done=true

    // Retry sequence (linear pipeline)
    subgraph->edges(retryDensity, retryCCDensity);
    subgraph->edges(retryCCDensity, retryVelocityFlux);
    subgraph->edges(retryVelocityFlux, retryHvac);
    subgraph->edges(retryHvac, retryInitDiv);
    subgraph->edges(retryInitDiv, retryDivPart1);
    subgraph->edges(retryDivPart1, retryDivExchange);
    subgraph->edges(retryDivExchange, retryDivPart2);
    subgraph->edges(retryDivPart2, retryPressure);
    subgraph->edges(retryPressure, retryVelocityPredictor);

    // End of retry sequence → RetryLoopState (check if another retry needed)
    subgraph->edges(retryVelocityPredictor, retryLoopSM);

    // RetryLoopState → cycle back to retry sequence OR exit
    // Both edges exist, but only one path will be taken based on done flag:
    // - done=false: cycle to retryDensity (which processes and loops back)
    // - done=true: passthrough to retryExit (which emits MeshData)
    subgraph->edges(retryLoopSM, retryDensity);  // Cycle back for another retry
    subgraph->edges(retryLoopSM, retryExit);     // Exit path

    // Graph output: RetryExitTask converts RetrySequenceData to MeshData
    subgraph->outputs(retryExit);

    return subgraph;
}

#endif // CHANGE_TIMESTEP_SUBGRAPH_H
