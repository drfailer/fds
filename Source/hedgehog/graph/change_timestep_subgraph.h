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
/// Dataflow:
///   CheckRetryTask → RetryDensity → ... → RetryVelocityPredictor → RetryLoopSM
///
/// RetryLoopSM uses type-based routing for two output types:
///   - RetrySequenceData → cycles back to RetryDensity for another retry
///   - MeshData → exits the sub-graph (retry complete or no retry needed)
///
/// When no retry is needed, CheckRetryTask sets done=true and all pipeline
/// tasks pass the data through without processing. RetryLoopState then emits
/// MeshData tokens to exit.
///
/// Termination is data-driven: RetryLoopState records (t, dt) from each token
/// it processes. canTerminate() returns true when t + dt >= tEnd, which keeps
/// the cycle alive across all time steps until the simulation ends.
///
/// @param tEnd Simulation end time (used for data-driven cycle termination)
/// @return Shared pointer to the constructed sub-graph
inline auto buildChangeTimeStepSubgraph(double tEnd) {
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

    // --- Create retry loop state manager (data-driven canTerminate) ---
    auto retryLoopSM = std::make_shared<RetryLoopStateManager>(
        std::make_shared<RetryLoopState>(tEnd), "RetryLoop");

    // --- Wire the sub-graph ---

    // Entry point: CheckRetryTask receives BarrierData
    subgraph->inputs(checkRetry);

    // CheckRetryTask always sends into the pipeline (done=true data passes through)
    subgraph->edges(checkRetry, retryDensity);

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

    // End of retry sequence → RetryLoopState
    subgraph->edges(retryVelocityPredictor, retryLoopSM);

    // Cycle: RetryLoopState emits RetrySequenceData → back to retryDensity
    subgraph->edges(retryLoopSM, retryDensity);

    // Exit: RetryLoopState emits MeshData → subgraph output
    subgraph->outputs(retryLoopSM);

    return subgraph;
}

#endif // CHANGE_TIMESTEP_SUBGRAPH_H
