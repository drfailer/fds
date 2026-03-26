#ifndef PRESSURE_ITERATION_SUBGRAPH_H
#define PRESSURE_ITERATION_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/pressure_iteration_data.h"
#include "../data/mesh_data.h"
#include "../data/termination_data.h"
#include "../task/baroclinic_kernel_task.h"
#include "../task/pressure_iteration_tasks.h"
#include "../task/mesh_exchange_task.h"
#include "../state/pressure_convergence_state.h"
#include "../state/barrier_state.h"
#include "../tool/mesh_dependency_graph.h"
#include "mesh_exchange_graph.h"

/// Build the pressure iteration sub-graph.
///
/// Architecture (parallel pull-only exchange via MeshExchangeGraph):
///
///   BaroclinicKernel(task, parallel)
///     -> MeshExchangeGraph (dependency-managed parallel exchange cycle)
///     -> PressureSolveKernel(task, parallel)
///     -> PressureConvergence(barrier: exchange(5) + vel_error + check)
///       -> PressureIterMeshData -> BaroclinicKernel (outer cycle)
///       -> MeshData -> subgraph output (converged)
///
/// IMPORTANT: fds_pressure_iteration_init() and the first
/// fds_pressure_iteration_increment() must be called in the upstream
/// barrier BEFORE entering this subgraph.
///
/// @param nmeshes Number of local meshes
/// @param kernelThreads Number of threads for parallel kernel tasks
/// @param exchangeThreads Number of threads for the flux exchange task
/// @param predictor True for predictor phase
/// @param depGraph Pre-built mesh dependency graph
/// @param commService Pointer to the MPI comm service (unused, kept for API compat)
/// @param presFlag Pressure solver flag (FFT_FLAG=0, ULMAT_FLAG=3)
inline auto buildPressureIterationSubgraph(int nmeshes,
                                            size_t kernelThreads,
                                            size_t exchangeThreads,
                                            bool predictor,
                                            std::shared_ptr<MeshDependencyGraph> depGraph,
                                            [[maybe_unused]] void *commService,
                                            int presFlag = 0) {
    using SubGraphType = hh::Graph<2, MeshData, TerminationData, MeshData>;
    auto subgraph = std::make_shared<SubGraphType>("PressureIteration");

    // --- Parallel kernel tasks ---
    auto baroclinicKernel = std::make_shared<BaroclinicKernelTask>(kernelThreads);
    auto solveKernel = std::make_shared<PressureSolveKernelTask>(kernelThreads, presFlag);

    // --- Dependency-aware parallel exchange (encapsulated cycle) ---
    auto exchangeGraph = std::make_shared<MeshExchangeGraph<MeshData>>(
        depGraph, std::make_shared<FluxExchangeTask>(exchangeThreads));

    // --- Convergence barrier (exchange2 + vel_error + convergence check) ---
    auto convergenceSM = std::make_shared<PressureConvergenceManager>(
        std::make_shared<PressureConvergenceState>(nmeshes, predictor),
        "PressureConvergence");

    // --- Wire the sub-graph ---

    // Entry: MeshData -> baroclinicKernel
    subgraph->input<MeshData>(baroclinicKernel);

    // TerminationData -> exchange graph + convergence
    MeshExchangeGraph<MeshData>::wireTermination(subgraph, exchangeGraph);
    subgraph->input<TerminationData>(convergenceSM);

    // baroclinicKernel -> MeshExchangeGraph -> PressureSolveKernel
    subgraph->edges(baroclinicKernel, exchangeGraph);
    subgraph->edges(exchangeGraph, solveKernel);

    // PressureSolveKernel -> PressureConvergence (barrier)
    subgraph->edges(solveKernel, convergenceSM);

    // Outer cycle: PressureIterMeshData -> back to baroclinicKernel
    subgraph->edge<PressureIterMeshData>(convergenceSM, baroclinicKernel);

    // Exit: MeshData -> subgraph output
    subgraph->outputs(convergenceSM);

    return subgraph;
}

#endif // PRESSURE_ITERATION_SUBGRAPH_H
