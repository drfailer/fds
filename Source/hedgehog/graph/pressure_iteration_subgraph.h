#ifndef PRESSURE_ITERATION_SUBGRAPH_H
#define PRESSURE_ITERATION_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/termination_data.h"
#include "../task/pressure_parallel_task.h"
#include "../state/pressure_convergence_state.h"
#include "../state/post_exchange_router_state.h"
#include "../state/barrier_state.h"
#include "../tool/mesh_dependency_graph.h"
#include "../tool/thread_budget.h"
#include "mesh_exchange_graph.h"

/// Build the pressure iteration sub-graph.
///
/// Shared between predictor and corrector phases. Phase routing is handled
/// internally via MeshData::phase (0=predictor, 1=corrector).
///
/// Inputs:
///   - MeshData<PredictorPressure>: from predictor pipeline
///   - MeshData<CorrectorPressure>: from corrector pipeline
///   - TerminationData: shutdown signal
///
/// Outputs:
///   - MeshData<PredictorPressure>: converged predictor result
///   - MeshData<CorrectorPressure>: converged corrector result
///
/// The 3 per-mesh parallel kernels (Baroclinic, Solve, VelocityError) are
/// packed into a single PressureParallelTask thread pool. Between phases,
/// exchange operations route data via distinct MeshState tags:
///   Baroclinic → PreSolveExch → exchange → SolvePhase → Solve
///   Solve → PostSolveExch → exchange → VelErrorPhase → VelError
///   VelError → Pressure → Convergence → cycle or exit
///
/// Two wiring modes:
///
/// **Single-process mode** (commService == nullptr):
///   PressureParallel → ExchangeInputRetag → MeshExchangeGraph<Pressure>
///     → PostExchangeRouter → PressureParallel (SolvePhase/VelErrorPhase)
///
/// **MPI mode** (commService != nullptr):
///   PressureParallel → RetaggingBarrier<PreSolveExch,SolvePhase> → PressureParallel
///   PressureParallel → RetaggingBarrier<PostSolveExch,VelErrorPhase> → PressureParallel
///
/// IMPORTANT: fds_pressure_iteration_init() and the first
/// fds_pressure_iteration_increment() must be called in the upstream
/// barrier BEFORE entering this subgraph.
///
/// @param nmeshes Number of local meshes
/// @param budget Thread budget for task thread allocation
/// @param depGraph Pre-built mesh dependency graph
/// @param commService Pointer to the MPI comm service (passed to MeshExchangeGraph)
/// @param presFlag Pressure solver flag (FFT_FLAG=0, ULMAT_FLAG=3)
inline auto buildPressureIterationSubgraph(int nmeshes,
                                            const ThreadBudget &budget,
                                            std::shared_ptr<MeshDependencyGraph> depGraph,
                                            hh::comm::CommService *commService,
                                            int presFlag = 0) {
    using SubGraphType = hh::Graph<3,
        MeshData<MeshState::PredictorPressure>,
        MeshData<MeshState::CorrectorPressure>,
        TerminationData,
        MeshData<MeshState::PredictorPressure>,
        MeshData<MeshState::CorrectorPressure>>;
    auto subgraph = std::make_shared<SubGraphType>("PressureIteration");

    // --- Packed parallel task (threads from budget) ---
    auto pressureParallelTask = std::make_shared<PressureParallelTask>(
        budget.pressureParallel, presFlag);

    // --- Convergence barrier (convergence check only) ---
    auto convergenceSM = std::make_shared<PressureConvergenceManager>(
        std::make_shared<PressureConvergenceState>(nmeshes),
        "PressureConvergence");

    // Entry: PredPressure + CorrPressure -> pressureParallelTask (baroclinic phase)
    subgraph->template input<MeshData<MeshState::PredictorPressure>>(pressureParallelTask);
    subgraph->template input<MeshData<MeshState::CorrectorPressure>>(pressureParallelTask);
    subgraph->template input<TerminationData>(convergenceSM);
    subgraph->template input<TerminationData>(pressureParallelTask);

    if (commService) {
        // --- MPI mode: retagging barriers for exchange ---
        auto preSolveExchange = makeRetaggingBarrier<MeshState::PreSolveExch, MeshState::SolvePhase>(
            nmeshes, "PreSolveExchange",
            "MESH_EXCHANGE(5)",
            [](auto&) { fds_mesh_exchange(5); });

        auto postSolveExchange = makeRetaggingBarrier<MeshState::PostSolveExch, MeshState::VelErrorPhase>(
            nmeshes, "PostSolveExchange",
            "MESH_EXCHANGE(5)",
            [](auto&) { fds_mesh_exchange(5); });

        // PressureParallel ↔ PreSolveExchange ↔ PressureParallel
        subgraph->edges(pressureParallelTask, preSolveExchange);
        subgraph->edges(preSolveExchange, pressureParallelTask);

        // PressureParallel ↔ PostSolveExchange ↔ PressureParallel
        subgraph->edges(pressureParallelTask, postSolveExchange);
        subgraph->edges(postSolveExchange, pressureParallelTask);

    } else {
        // --- Single-process mode: retag → exchange graph → router ---
        auto retagTask = std::make_shared<ExchangeInputRetagTask>();

        auto exchange = std::make_shared<MeshExchangeGraph<MeshState::Pressure>>(
            depGraph, budget.exchangePush, budget.exchangePull,
            std::vector<int>{5},
            nullptr, "Exchange");

        auto routerSM = std::make_shared<PostExchangeRouterManager>(
            std::make_shared<PostExchangeRouterState>(),
            "PostExchangeRouter");

        // TerminationData -> router (for canTerminate to break inner cycle)
        subgraph->template input<TerminationData>(routerSM);

        // PressureParallel → retag(PreSolveExch+PostSolveExch → Pressure)
        //   → Exchange → Router(SolvePhase+VelErrorPhase) → PressureParallel
        subgraph->edges(pressureParallelTask, retagTask);
        subgraph->edges(retagTask, exchange);
        subgraph->edges(exchange, routerSM);
        subgraph->edges(routerSM, pressureParallelTask);
    }

    // VelError output (Pressure) → Convergence
    subgraph->edges(pressureParallelTask, convergenceSM);

    // Cycle: MeshData<Pressure> → back to pressureParallelTask (baroclinic phase)
    // Typed edge: avoid routing PredPressure/CorrPressure back into cycle
    subgraph->template edge<MeshData<MeshState::Pressure>>(convergenceSM, pressureParallelTask);

    // Exit: PredPressure + CorrPressure → subgraph output
    subgraph->outputs(convergenceSM);

    return subgraph;
}

#endif // PRESSURE_ITERATION_SUBGRAPH_H
