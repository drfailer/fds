#ifndef PRESSURE_ITERATION_SUBGRAPH_H
#define PRESSURE_ITERATION_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/pressure_iteration_data.h"
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

/// Lightweight wrapper task: MeshData -> SolvePhaseData.
/// Used in MPI mode (barrier-based exchange) and single-process mode
/// (PostExchangeRouter) to route data after exchange.
class SolvePhaseWrapperTask
    : public hh::AbstractTask<1, MeshData, SolvePhaseData> {
public:
    SolvePhaseWrapperTask()
        : hh::AbstractTask<1, MeshData, SolvePhaseData>(
              "SolvePhaseWrap", 1) {}
    void execute(std::shared_ptr<MeshData> md) override {
        this->addResult(std::make_shared<SolvePhaseData>(std::move(md)));
    }
    std::shared_ptr<hh::AbstractTask<1, MeshData, SolvePhaseData>>
    copy() override { return std::make_shared<SolvePhaseWrapperTask>(); }
};

/// Lightweight wrapper task: MeshData -> VelErrorPhaseData.
class VelErrorPhaseWrapperTask
    : public hh::AbstractTask<1, MeshData, VelErrorPhaseData> {
public:
    VelErrorPhaseWrapperTask()
        : hh::AbstractTask<1, MeshData, VelErrorPhaseData>(
              "VelErrorPhaseWrap", 1) {}
    void execute(std::shared_ptr<MeshData> md) override {
        this->addResult(std::make_shared<VelErrorPhaseData>(std::move(md)));
    }
    std::shared_ptr<hh::AbstractTask<1, MeshData, VelErrorPhaseData>>
    copy() override { return std::make_shared<VelErrorPhaseWrapperTask>(); }
};

/// Build the pressure iteration sub-graph.
///
/// Two wiring modes:
///
/// **Single-process mode** (commService == nullptr):
/// Single exchange graph, double-buffered, pipeline-parallel with cycle.
///
///   BaroclinicKernel(task, sets exchangeRound=0)
///     -> MeshExchangeGraph (round 0: pre-solve)
///     -> PostExchangeRouter(SM, pass-through)
///       -> SolvePhaseData -> PressureSolveKernel(task, sets exchangeRound=1)
///                              -> MeshExchangeGraph (round 1: post-solve, SAME graph)
///                              -> PostExchangeRouter
///       -> VelErrorPhaseData -> VelocityErrorTask(task)
///                                 -> PressureConvergence(barrier)
///                                   -> PressureIterMeshData -> BaroclinicKernel (cycle)
///                                   -> MeshData -> subgraph output (converged)
///
/// **MPI mode** (commService != nullptr):
/// Barrier-based exchange using fds_mesh_exchange(5). Avoids spawning
/// CommunicatorTask daemon threads inside the pressure iteration, which
/// cause severe CPU spinning and performance degradation.
///
///   BaroclinicKernel -> PreExchangeBarrier(fds_mesh_exchange(5))
///     -> [wrap] -> SolveKernel -> PostExchangeBarrier(fds_mesh_exchange(5))
///     -> [wrap] -> VelocityError -> Convergence
///       -> PressureIterMeshData -> BaroclinicKernel (cycle)
///       -> MeshData -> subgraph output (converged)
///
/// IMPORTANT: fds_pressure_iteration_init() and the first
/// fds_pressure_iteration_increment() must be called in the upstream
/// barrier BEFORE entering this subgraph.
///
/// @param nmeshes Number of local meshes
/// @param budget Thread budget for task thread allocation
/// @param predictor True for predictor phase
/// @param depGraph Pre-built mesh dependency graph
/// @param commService Pointer to the MPI comm service (passed to MeshExchangeGraph)
/// @param presFlag Pressure solver flag (FFT_FLAG=0, ULMAT_FLAG=3)
inline auto buildPressureIterationSubgraph(int nmeshes,
                                            const ThreadBudget &budget,
                                            bool predictor,
                                            std::shared_ptr<MeshDependencyGraph> depGraph,
                                            hh::comm::CommService *commService,
                                            int presFlag = 0) {
    using SubGraphType = hh::Graph<2, MeshData, TerminationData, MeshData>;
    auto subgraph = std::make_shared<SubGraphType>("PressureIteration");

    // --- Parallel kernel tasks (threads from budget) ---
    auto baroclinicKernel = std::make_shared<BaroclinicKernelTask>(budget.baroclinic);
    auto solveKernel = std::make_shared<PressureSolveKernelTask>(budget.pressureSolve, presFlag);
    auto velErrorTask = std::make_shared<VelocityErrorTask>(budget.velError);

    // --- Convergence barrier (convergence check only) ---
    auto convergenceSM = std::make_shared<PressureConvergenceManager>(
        std::make_shared<PressureConvergenceState>(nmeshes, predictor),
        "PressureConvergence");

    // Entry: MeshData -> baroclinicKernel
    subgraph->input<MeshData>(baroclinicKernel);
    subgraph->input<TerminationData>(convergenceSM);

    if (commService) {
        // --- MPI mode: barrier-based exchange, linear pipeline ---
        // Uses fds_mesh_exchange(5) barriers instead of CommunicatorTask.
        // CommunicatorTask daemon threads busy-wait on MPI_Iprobe in a tight
        // loop, consuming 100% CPU and starving compute threads. Barrier-based
        // exchange avoids this overhead entirely.

        auto preSolveExchange = makeBarrierSM(nmeshes, "PreSolveExchange",
            "MESH_EXCHANGE(5)", [](std::vector<std::shared_ptr<MeshData>> &) {
                fds_mesh_exchange(5);
            });

        auto postSolveExchange = makeBarrierSM(nmeshes, "PostSolveExchange",
            "MESH_EXCHANGE(5)", [](std::vector<std::shared_ptr<MeshData>> &) {
                fds_mesh_exchange(5);
            });

        auto solveWrap = std::make_shared<SolvePhaseWrapperTask>();
        auto velErrorWrap = std::make_shared<VelErrorPhaseWrapperTask>();

        // BaroclinicKernel -> PreSolveExchange barrier -> wrap -> SolveKernel
        subgraph->edges(baroclinicKernel, preSolveExchange);
        subgraph->edges(preSolveExchange, solveWrap);
        subgraph->edges(solveWrap, solveKernel);

        // SolveKernel -> PostSolveExchange barrier -> wrap -> VelocityError
        subgraph->edges(solveKernel, postSolveExchange);
        subgraph->edges(postSolveExchange, velErrorWrap);
        subgraph->edges(velErrorWrap, velErrorTask);

        // VelErrorTask -> Convergence
        subgraph->edges(velErrorTask, convergenceSM);

    } else {
        // --- Single-process mode: single exchange graph with cycle ---
        // Pipeline-parallel via double-buffered state.

        auto exchange = std::make_shared<MeshExchangeGraph>(
            depGraph, budget.exchangePush, budget.exchangePull,
            std::vector<int>{5},
            nullptr, "Exchange");

        auto routerSM = std::make_shared<PostExchangeRouterManager>(
            std::make_shared<PostExchangeRouterState>(),
            "PostExchangeRouter");

        // TerminationData -> router (for canTerminate to break inner cycle)
        subgraph->input<TerminationData>(routerSM);

        // BaroclinicKernel -> Exchange (round 0)
        subgraph->edges(baroclinicKernel, exchange);

        // Exchange -> Router
        subgraph->edges(exchange, routerSM);

        // Router dispatches by exchangeRound:
        //   round 0 -> SolvePhaseData -> solveKernel
        //   round 1 -> VelErrorPhaseData -> velErrorTask
        subgraph->edge<SolvePhaseData>(routerSM, solveKernel);
        subgraph->edge<VelErrorPhaseData>(routerSM, velErrorTask);

        // SolveKernel -> Exchange (round 1, same graph)
        subgraph->edges(solveKernel, exchange);

        // VelErrorTask -> Convergence
        subgraph->edges(velErrorTask, convergenceSM);
    }

    // Outer cycle: PressureIterMeshData -> back to baroclinicKernel
    subgraph->edge<PressureIterMeshData>(convergenceSM, baroclinicKernel);

    // Exit: MeshData -> subgraph output
    subgraph->outputs(convergenceSM);

    return subgraph;
}

#endif // PRESSURE_ITERATION_SUBGRAPH_H
