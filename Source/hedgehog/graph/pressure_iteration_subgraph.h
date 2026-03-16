#ifndef PRESSURE_ITERATION_SUBGRAPH_H
#define PRESSURE_ITERATION_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/pressure_iteration_data.h"
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/termination_signal.h"
#include "../task/pressure_iteration_tasks.h"
#include "../state/pressure_iteration_state.h"

/// Build the pressure iteration sub-graph.
///
/// 4-node architecture:
///   PressurePreKernel -> PressureSolveKernel -> PressureSolveCollector -> PressurePostLoopSM
///                      ^                                                   | (cycle)
///                      +---------------------------------------------------+
///
/// PressurePreKernelTask accepts BarrierData (first entry) and
/// PressureIterData (cycle). Runs Phase 1 (baroclinic + exchange),
/// then scatters MeshData for parallel Phase 2 kernel.
///
/// PressureSolveCollector gathers N MeshData into PressureIterData.
///
/// PressurePostLoopSM runs Phase 3 (MESH_EXCHANGE(5) + velocity error +
/// convergence check), then uses type-based routing:
///   - PressureIterData -> cycles back to PressurePreKernel
///   - MeshData -> exits the sub-graph (converged)
///
/// Termination: The cycle stays alive across all time steps. canTerminate()
/// uses reachedEnd() && lastConverged() (data-driven, fires after final
/// convergence on last timestep) with TerminationSignal as fallback.
///
/// Only FFT solver is supported in parallel mode. ULMAT/GLMAT/UGLMAT
/// and CC_IBM cases fall back to sequential PressureIterationTask.
///
/// @param tEnd Simulation end time
/// @param nmeshes Number of local meshes
/// @param kernelThreads Number of threads for parallel kernel tasks
/// @param predictor True for predictor phase (calls init_change_time_step on exit)
/// @param termSignal Shared termination signal from the main timestep loop
/// @return Shared pointer to the constructed sub-graph
inline auto buildPressureIterationSubgraph(double tEnd, int nmeshes,
                                            size_t kernelThreads,
                                            bool predictor,
                                            std::shared_ptr<TerminationSignal> termSignal) {
    using SubGraphType = hh::Graph<1, BarrierData, MeshData>;
    auto subgraph = std::make_shared<SubGraphType>("PressureIteration");

    auto preKernel = std::make_shared<PressurePreKernelTask>();
    auto solveKernel = std::make_shared<PressureSolveKernelTask>(kernelThreads);
    auto solveCollSM = std::make_shared<hh::StateManager<
        1, MeshData, PressureIterData>>(
        std::make_shared<PressureSolveCollector>(nmeshes), "PressureSolveCollector");
    auto postLoopSM = std::make_shared<PressurePostLoopStateManager>(
        std::make_shared<PressurePostLoopState>(tEnd, predictor, termSignal),
        "PressurePostLoop");

    // Entry: BarrierData -> PreKernel
    subgraph->inputs(preKernel);

    // PreKernel -> parallel kernel (MeshData scatter)
    subgraph->edges(preKernel, solveKernel);

    // Kernel -> collector (MeshData -> PressureIterData)
    subgraph->edges(solveKernel, solveCollSM);

    // Collector -> post-loop (PressureIterData)
    subgraph->edges(solveCollSM, postLoopSM);

    // Cycle: PressureIterData -> back to PreKernel
    subgraph->edges(postLoopSM, preKernel);

    // Exit: MeshData -> subgraph output
    subgraph->outputs(postLoopSM);

    return subgraph;
}

#endif // PRESSURE_ITERATION_SUBGRAPH_H
