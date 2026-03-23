#ifndef PRESSURE_ITERATION_SUBGRAPH_H
#define PRESSURE_ITERATION_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/pressure_iteration_data.h"
#include "../data/mesh_data.h"
#include "../data/termination_signal.h"
#include "../task/pressure_iteration_tasks.h"
#include "../state/pressure_iteration_state.h"

/// Build the pressure iteration sub-graph.
///
/// 2-state + 1-task architecture:
///   PreCollector(state) -> SolveKernel(task) -> PostLoop(state)
///                ^                                  | (cycle)
///                +----------------------------------+
///
/// PressurePreCollector accepts MeshData (initial entry, collects N tokens)
/// and PressureIterData (cycle). Runs Phase 1 (pressure iteration init +
/// baroclinic correction + mesh exchange), then scatters MeshData for
/// parallel Phase 2 kernel.
///
/// PressureSolveKernelTask runs Phase 2 per mesh in parallel.
///
/// PressurePostLoop collects N MeshData from the solve kernel,
/// runs Phase 3 (exchange + velocity error + convergence), then routes
/// by type (PressureIterData -> cycle, MeshData -> exit).
///
/// Type-specific cycle edge (edge<PressureIterData>) prevents exit MeshData
/// from cycling back to the PreCollector.
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
    using SubGraphType = hh::Graph<1, MeshData, MeshData>;
    auto subgraph = std::make_shared<SubGraphType>("PressureIteration");

    auto preCollSM = std::make_shared<PressurePreCollectorManager>(
        std::make_shared<PressurePreCollector>(nmeshes),
        "PressurePreCollector");
    auto solveKernel = std::make_shared<PressureSolveKernelTask>(kernelThreads, presFlag);
    auto postLoopSM = std::make_shared<PressurePostLoopManager>(
        std::make_shared<PressurePostLoop>(nmeshes, tEnd, predictor, termSignal),
        "PressurePostLoop");

    // Entry: MeshData -> PreCollector (collects N tokens, then Phase 1)
    subgraph->inputs(preCollSM);

    // PreCollector -> parallel kernel (MeshData scatter)
    subgraph->edges(preCollSM, solveKernel);

    // Kernel -> PostLoop (MeshData collect N -> Phase 3 -> route)
    subgraph->edges(solveKernel, postLoopSM);

    // Cycle: PressureIterData ONLY -> back to PreCollector.
    // Uses edge<> (single-type) instead of edges() (all-types) to prevent
    // PostLoop's exit MeshData from cycling back to PreCollector's MeshData input.
    subgraph->edge<PressureIterData>(postLoopSM, preCollSM);

    // Exit: MeshData -> subgraph output
    subgraph->outputs(postLoopSM);

    return subgraph;
}

#endif // PRESSURE_ITERATION_SUBGRAPH_H
