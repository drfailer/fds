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
/// 3-node architecture (was 4 — SolveCollector merged into PostCollector):
///   PressurePreKernel -> PressureSolveKernel -> PressurePostCollector
///                      ^                                | (cycle)
///                      +--------------------------------+
///
/// PressurePreKernelTask accepts BarrierData (first entry) and
/// PressureIterData (cycle). Runs Phase 1 (baroclinic correction +
/// mesh exchange), then scatters MeshData for parallel Phase 2 kernel.
///
/// PressureSolveKernelTask runs Phase 2 per mesh in parallel, including
/// match_velocity_flux_kernel (moved from PreKernel for parallelization).
///
/// PressurePostCollector merges the former SolveCollector + PostLoopState:
/// collects N MeshData, runs Phase 3 (exchange + velocity error + convergence),
/// then routes by type (PressureIterData -> cycle, MeshData -> exit).
///
/// @param tEnd Simulation end time
/// @param nmeshes Number of local meshes
/// @param kernelThreads Number of threads for parallel kernel tasks
/// @param predictor True for predictor phase (calls init_change_time_step on exit)
/// @param termSignal Shared termination signal from the main timestep loop
/// @param presFlag Pressure solver flag (FFT_FLAG=0, ULMAT_FLAG=3)
/// @return Shared pointer to the constructed sub-graph
inline auto buildPressureIterationSubgraph(double tEnd, int nmeshes,
                                            size_t kernelThreads,
                                            bool predictor,
                                            std::shared_ptr<TerminationSignal> termSignal,
                                            int presFlag = 0) {
    using SubGraphType = hh::Graph<1, BarrierData, MeshData>;
    auto subgraph = std::make_shared<SubGraphType>("PressureIteration");

    auto preKernel = std::make_shared<PressurePreKernelTask>();
    auto solveKernel = std::make_shared<PressureSolveKernelTask>(kernelThreads, presFlag);
    auto postCollSM = std::make_shared<PressurePostCollectorManager>(
        std::make_shared<PressurePostCollector>(nmeshes, tEnd, predictor, termSignal),
        "PressurePostCollector");

    // Entry: BarrierData -> PreKernel
    subgraph->inputs(preKernel);

    // PreKernel -> parallel kernel (MeshData scatter)
    subgraph->edges(preKernel, solveKernel);

    // Kernel -> merged collector+post-loop (MeshData -> collect N -> route)
    subgraph->edges(solveKernel, postCollSM);

    // Cycle: PressureIterData -> back to PreKernel
    subgraph->edges(postCollSM, preKernel);

    // Exit: MeshData -> subgraph output
    subgraph->outputs(postCollSM);

    return subgraph;
}

#endif // PRESSURE_ITERATION_SUBGRAPH_H
