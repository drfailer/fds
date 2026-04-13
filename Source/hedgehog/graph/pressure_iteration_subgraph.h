#ifndef PRESSURE_ITERATION_SUBGRAPH_H
#define PRESSURE_ITERATION_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/termination_data.h"
#include "../task/baroclinic_kernel_task.h"
#include "../task/pressure_iteration_tasks.h"
#include "../task/velocity_error_task.h"
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
/// Two wiring modes:
///
/// **Single-process mode** (commService == nullptr):
/// Single exchange graph, double-buffered, pipeline-parallel with cycle.
///
///   BaroclinicKernel(task, 3 inputs: Pred+Corr+Pressure)
///     -> MeshExchangeGraph<Pressure> (round 0: pre-solve)
///     -> PostExchangeRouter(SM, pass-through)
///       -> MeshData<SolvePhase> -> PressureSolveKernel(task, sets exchangeRound=1)
///                                    -> MeshExchangeGraph (round 1, SAME graph)
///                                    -> PostExchangeRouter
///       -> MeshData<VelErrorPhase> -> VelocityErrorTask(task)
///                                       -> PressureConvergence(barrier)
///                                         -> MeshData<Pressure> -> BaroclinicKernel (cycle)
///                                         -> MeshData<Pred|CorrPressure> -> output
///
/// **MPI mode** (commService != nullptr):
/// Barrier-based exchange using fds_mesh_exchange(5).
///
///   BaroclinicKernel -> PreExchangeBarrier<Pressure>(fds_mesh_exchange(5))
///     -> SolveKernel -> PostExchangeBarrier<Pressure>(fds_mesh_exchange(5))
///     -> VelocityError -> Convergence
///       -> MeshData<Pressure> -> BaroclinicKernel (cycle)
///       -> MeshData<Pred|CorrPressure> -> output
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

    // --- Parallel kernel tasks (threads from budget) ---
    auto baroclinicKernel = std::make_shared<BaroclinicKernelTask>(budget.baroclinic);
    auto solveKernel = std::make_shared<PressureSolveKernelTask>(budget.pressureSolve, presFlag);
    auto velErrorTask = std::make_shared<VelocityErrorTask>(budget.velError);

    // --- Convergence barrier (convergence check only) ---
    auto convergenceSM = std::make_shared<PressureConvergenceManager>(
        std::make_shared<PressureConvergenceState>(nmeshes),
        "PressureConvergence");

    // Entry: PredPressure + CorrPressure -> baroclinicKernel
    subgraph->template input<MeshData<MeshState::PredictorPressure>>(baroclinicKernel);
    subgraph->template input<MeshData<MeshState::CorrectorPressure>>(baroclinicKernel);
    subgraph->template input<TerminationData>(convergenceSM);

    if (commService) {
        // --- MPI mode: barrier-based exchange, linear pipeline ---
        // Uses fds_mesh_exchange(5) barriers instead of CommunicatorTask.
        // CommunicatorTask daemon threads busy-wait on MPI_Iprobe in a tight
        // loop, consuming 100% CPU and starving compute threads. Barrier-based
        // exchange avoids this overhead entirely.

        auto preSolveExchange = makeBarrierSM<MeshState::Pressure>(
            nmeshes, "PreSolveExchange",
            "MESH_EXCHANGE(5)",
            [](std::vector<std::shared_ptr<MeshData<MeshState::Pressure>>> &) {
                fds_mesh_exchange(5);
            });

        auto postSolveExchange = makeBarrierSM<MeshState::Pressure>(
            nmeshes, "PostSolveExchange",
            "MESH_EXCHANGE(5)",
            [](std::vector<std::shared_ptr<MeshData<MeshState::Pressure>>> &) {
                fds_mesh_exchange(5);
            });

        // BaroclinicKernel -> PreSolveExchange barrier -> SolveKernel
        subgraph->edges(baroclinicKernel, preSolveExchange);
        subgraph->edges(preSolveExchange, solveKernel);

        // SolveKernel -> PostSolveExchange barrier -> VelocityError
        subgraph->edges(solveKernel, postSolveExchange);
        subgraph->edges(postSolveExchange, velErrorTask);

        // VelErrorTask -> Convergence
        subgraph->edges(velErrorTask, convergenceSM);

    } else {
        // --- Single-process mode: single exchange graph with cycle ---
        // Pipeline-parallel via double-buffered state.

        auto exchange = std::make_shared<MeshExchangeGraph<MeshState::Pressure>>(
            depGraph, budget.exchangePush, budget.exchangePull,
            std::vector<int>{5},
            nullptr, "Exchange");

        auto routerSM = std::make_shared<PostExchangeRouterManager>(
            std::make_shared<PostExchangeRouterState>(),
            "PostExchangeRouter");

        // TerminationData -> router (for canTerminate to break inner cycle)
        subgraph->template input<TerminationData>(routerSM);

        // BaroclinicKernel -> Exchange (round 0)
        subgraph->edges(baroclinicKernel, exchange);

        // Exchange -> Router
        subgraph->edges(exchange, routerSM);

        // Router dispatches by exchangeRound:
        //   round 0 -> MeshData<SolvePhase> -> solveKernel
        //   round 1 -> MeshData<VelErrorPhase> -> velErrorTask
        subgraph->template edge<MeshData<MeshState::SolvePhase>>(routerSM, solveKernel);
        subgraph->template edge<MeshData<MeshState::VelErrorPhase>>(routerSM, velErrorTask);

        // SolveKernel -> Exchange (round 1, same graph)
        subgraph->edges(solveKernel, exchange);

        // VelErrorTask -> Convergence
        subgraph->edges(velErrorTask, convergenceSM);
    }

    // Cycle: MeshData<Pressure> -> back to baroclinicKernel
    subgraph->template edge<MeshData<MeshState::Pressure>>(convergenceSM, baroclinicKernel);

    // Exit: PredPressure + CorrPressure -> subgraph output
    subgraph->outputs(convergenceSM);

    return subgraph;
}

#endif // PRESSURE_ITERATION_SUBGRAPH_H
