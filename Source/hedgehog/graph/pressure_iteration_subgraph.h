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
#include "../state/barrier_state.h"
#include "../tool/mesh_dependency_graph.h"
#include "../tool/thread_budget.h"
#include "mesh_exchange_graph.h"

/// Build the pressure iteration sub-graph.
///
/// Architecture (push/buffer/pull exchange via MeshExchangeGraph):
///
///   BaroclinicKernel(task, parallel)
///     -> PreSolveExchange (push → deps gate → pull)
///     -> PressureSolveKernel(task, parallel)
///     -> PostSolveExchange (push → deps gate → pull)
///     -> VelocityErrorTask(task, parallel)
///     -> PressureConvergence(barrier: convergence check only)
///       -> PressureIterMeshData -> BaroclinicKernel (outer cycle)
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

    // --- Push/buffer/pull exchanges (linear pipelines, no internal cycle) ---
    auto preSolveExchange = std::make_shared<MeshExchangeGraph<MeshData>>(
        depGraph, budget.exchangePush, budget.exchangePull,
        commService, "PreSolveExchange");
    auto postSolveExchange = std::make_shared<MeshExchangeGraph<MeshData>>(
        depGraph, budget.exchangePush, budget.exchangePull,
        commService, "PostSolveExchange");

    // --- Convergence barrier (convergence check only) ---
    auto convergenceSM = std::make_shared<PressureConvergenceManager>(
        std::make_shared<PressureConvergenceState>(nmeshes, predictor),
        "PressureConvergence");

    // --- Wire the sub-graph ---

    // Entry: MeshData -> baroclinicKernel
    subgraph->input<MeshData>(baroclinicKernel);

    // TerminationData -> convergence (exchange graphs need none — no internal cycle)
    subgraph->input<TerminationData>(convergenceSM);

    // Pipeline: Baroclinic -> PreSolveExchange -> Solve -> PostSolveExchange -> VelError -> Convergence
    subgraph->edges(baroclinicKernel, preSolveExchange);
    subgraph->edges(preSolveExchange, solveKernel);
    subgraph->edges(solveKernel, postSolveExchange);
    subgraph->edges(postSolveExchange, velErrorTask);
    subgraph->edges(velErrorTask, convergenceSM);

    // Outer cycle: PressureIterMeshData -> back to baroclinicKernel
    subgraph->edge<PressureIterMeshData>(convergenceSM, baroclinicKernel);

    // Exit: MeshData -> subgraph output
    subgraph->outputs(convergenceSM);

    return subgraph;
}

#endif // PRESSURE_ITERATION_SUBGRAPH_H
