#ifndef PRESSURE_ITERATION_SUBGRAPH_H
#define PRESSURE_ITERATION_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/pressure_iteration_data.h"
#include "../data/mesh_data.h"
#include "../data/termination_signal.h"
#include "../task/baroclinic_kernel_task.h"
#include "../task/pressure_iteration_tasks.h"
#include "../state/pressure_convergence_state.h"
#include "../state/barrier_state.h"
#include "../state/exchange_orchestrator_state.h"
#include "../tool/mesh_dependency_graph.h"

/// Build the pressure iteration sub-graph.
///
/// Architecture (dependency-aware exchange + 2 parallel kernel tasks):
///
///   BaroclinicKernel(task, parallel)
///     → ExchangeOrchestrator(state: push copies + dependency gate)
///     → PressureSolveKernel(task, parallel)
///     → PressureConvergence(barrier: exchange(5) + vel_error + check)
///       → PressureIterMeshData → BaroclinicKernel (cycle)
///       → MeshData → subgraph output (converged)
///
/// The ExchangeOrchestrator replaces the global fds_mesh_exchange(5) barrier
/// with per-mesh dependency tracking.  When a mesh arrives, its data is
/// immediately pushed (copied) to all same-rank targets' OMESHes.  A mesh
/// can proceed to the pressure solve once all its receive-dependencies
/// have also pushed.
///
/// IMPORTANT: fds_pressure_iteration_init() and the first
/// fds_pressure_iteration_increment() must be called in the upstream
/// barrier BEFORE entering this subgraph.
///
/// @param tEnd Simulation end time
/// @param nmeshes Number of local meshes
/// @param kernelThreads Number of threads for parallel kernel tasks
/// @param predictor True for predictor phase
/// @param termSignal Shared termination signal
/// @param depGraph Pre-built mesh dependency graph
/// @param commService Pointer to the MPI comm service (unused, kept for API compat)
/// @param presFlag Pressure solver flag (FFT_FLAG=0, ULMAT_FLAG=3)
inline auto buildPressureIterationSubgraph(double tEnd, int nmeshes,
                                            size_t kernelThreads,
                                            bool predictor,
                                            std::shared_ptr<TerminationSignal> termSignal,
                                            std::shared_ptr<MeshDependencyGraph> depGraph,
                                            [[maybe_unused]] void *commService,
                                            int presFlag = 0) {
    using SubGraphType = hh::Graph<1, MeshData, MeshData>;
    auto subgraph = std::make_shared<SubGraphType>("PressureIteration");

    // --- Parallel kernel tasks ---
    auto baroclinicKernel = std::make_shared<BaroclinicKernelTask>(kernelThreads);
    auto solveKernel = std::make_shared<PressureSolveKernelTask>(kernelThreads, presFlag);

    // --- Dependency-aware exchange ---
    auto exchangeOrchestratorSM = std::make_shared<ExchangeOrchestratorManager>(
        std::make_shared<ExchangeOrchestratorState>(depGraph),
        "ExchangeOrchestrator");

    // --- Convergence barrier (exchange2 + vel_error + convergence check) ---
    auto convergenceSM = std::make_shared<PressureConvergenceManager>(
        std::make_shared<PressureConvergenceState>(
            nmeshes, tEnd, predictor, termSignal),
        "PressureConvergence");

    // --- Wire the sub-graph ---

    // Entry: MeshData -> baroclinicKernel
    subgraph->inputs(baroclinicKernel);

    // baroclinicKernel -> ExchangeOrchestrator (push copies + dependency gate)
    subgraph->edges(baroclinicKernel, exchangeOrchestratorSM);

    // ExchangeOrchestrator -> PressureSolveKernel (emits when deps satisfied)
    subgraph->edges(exchangeOrchestratorSM, solveKernel);

    // PressureSolveKernel -> PressureConvergence (barrier)
    subgraph->edges(solveKernel, convergenceSM);

    // Cycle: PressureIterMeshData -> back to baroclinicKernel
    subgraph->edge<PressureIterMeshData>(convergenceSM, baroclinicKernel);

    // Exit: MeshData -> subgraph output
    subgraph->outputs(convergenceSM);

    return subgraph;
}

#endif // PRESSURE_ITERATION_SUBGRAPH_H
