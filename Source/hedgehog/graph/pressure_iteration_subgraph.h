#ifndef PRESSURE_ITERATION_SUBGRAPH_H
#define PRESSURE_ITERATION_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/pressure_iteration_data.h"
#include "../data/mesh_data.h"
#include "../data/mesh_exchange_data.h"
#include "../data/termination_signal.h"
#include "../task/baroclinic_kernel_task.h"
#include "../task/pressure_iteration_tasks.h"
#include "../task/mesh_exchange_task.h"  // FluxExchangeTask
#include "../state/pressure_convergence_state.h"
#include "../state/barrier_state.h"
#include "../state/mesh_dependencies_manager_state.h"
#include "../tool/mesh_dependency_graph.h"

/// Build the pressure iteration sub-graph.
///
/// Architecture (dependency-aware bidirectional exchange):
///
///   BaroclinicKernel(task, parallel)
///     → MeshDependenciesManager(state: tracks arrivals, emits when deps met)
///     → MeshExchangeTask(task, single-threaded: bidirectional copies)
///     → PressureSolveKernel(task, parallel)
///     → PressureConvergence(barrier: exchange(5) + vel_error + check)
///       → PressureIterMeshData → BaroclinicKernel (cycle)
///       → MeshData → subgraph output (converged)
///
/// The MeshDependenciesManager state tracks which meshes have arrived from
/// the baroclinic kernel.  When a mesh's same-rank neighbors have all
/// arrived (or been exchanged), the state emits a MeshExchangeData token
/// containing the mesh and its list of non-yet-exchanged neighbors.
///
/// The MeshExchangeTask performs bidirectional copies (pull + push) for
/// each neighbor in the list, then emits the mesh as MeshData for the
/// pressure solve.
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

    // --- Dependency-aware bidirectional exchange ---
    auto depManagerSM = std::make_shared<MeshDependenciesManager>(
        std::make_shared<MeshDependenciesManagerState>(depGraph),
        "MeshDepsManager");
    auto exchangeTask = std::make_shared<FluxExchangeTask>();

    // --- Convergence barrier (exchange2 + vel_error + convergence check) ---
    auto convergenceSM = std::make_shared<PressureConvergenceManager>(
        std::make_shared<PressureConvergenceState>(
            nmeshes, tEnd, predictor, termSignal),
        "PressureConvergence");

    // --- Wire the sub-graph ---

    // Entry: MeshData -> baroclinicKernel
    subgraph->inputs(baroclinicKernel);

    // baroclinicKernel -> MeshDependenciesManager (arrival tracking)
    subgraph->edges(baroclinicKernel, depManagerSM);

    // MeshDependenciesManager -> MeshExchangeTask (bidirectional copies)
    subgraph->edges(depManagerSM, exchangeTask);

    // MeshExchangeTask -> PressureSolveKernel
    subgraph->edges(exchangeTask, solveKernel);

    // PressureSolveKernel -> PressureConvergence (barrier)
    subgraph->edges(solveKernel, convergenceSM);

    // Cycle: PressureIterMeshData -> back to baroclinicKernel
    subgraph->edge<PressureIterMeshData>(convergenceSM, baroclinicKernel);

    // Exit: MeshData -> subgraph output
    subgraph->outputs(convergenceSM);

    return subgraph;
}

#endif // PRESSURE_ITERATION_SUBGRAPH_H
