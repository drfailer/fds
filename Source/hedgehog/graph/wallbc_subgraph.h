#ifndef WALLBC_SUBGRAPH_H
#define WALLBC_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/wallbc_data.h"
#include "../task/wallbc_kernel_task.h"
#include "../state/wallbc_state.h"

/// Build the WallBC sub-graph (Pattern B: complex routine parallelization).
///
/// Three-phase architecture:
///   1. Sequential preprocessing (Orchestrator):
///      - Computes global parameters (DT_BC, CALL_HT_1D)
///      - WALL_BC_PREPROCESSING for each mesh:
///        * ASSIGN_GHOST_VALUE (OMESH reads for ghost cells)
///        * NEAR_SURFACE_GAS_VARIABLES_KERNEL (all wall cells)
///        * HEAT_TRANS_COEF (thermally-thick surfaces)
///
///   2. Parallel kernel execution (WallBCKernelTask):
///      - WALL_BC_PROCESS_CELLS_KERNEL (~90% of wall cells)
///      - Thread-safe routines: SURFACE_HEAT_TRANSFER, CALCULATE_ZZ_F, etc.
///      - Skips cells with cross-mesh dependencies
///
///   3. Sequential finalization (Collector):
///      - WALL_BC_FINALIZE for each mesh:
///        * HAS_BACK_MESH cells (thin walls spanning meshes)
///        * Thin wall lateral heat transfer
///        * Particle off-gassing via DEPOSIT_PARTICLE_MASS
///
/// This sub-graph replaces the sequential CorrWallBCTask.
///
/// @param nmeshes Number of meshes
/// @param kernelThreads Number of threads for parallel kernel task
/// @return Shared pointer to the constructed sub-graph
inline auto buildWallBCSubgraph(int nmeshes, size_t kernelThreads) {
    using SubGraphType = hh::Graph<1, MeshData, MeshData>;
    auto subgraph = std::make_shared<SubGraphType>("WallBC");

    // --- Create WallBC components ---
    auto wallBCOrchSM = std::make_shared<hh::StateManager<1, MeshData, WallBCWork>>(
        std::make_shared<WallBCOrchestrator>(nmeshes), "WallBCOrch");
    auto wallBCKernelTask = std::make_shared<WallBCKernelTask>(kernelThreads);
    auto wallBCCollectorSM = std::make_shared<hh::StateManager<1, WallBCWork, MeshData>>(
        std::make_shared<WallBCCollector>(nmeshes), "WallBCCollector");

    // --- Wire the sub-graph ---

    // Entry point: Orchestrator receives MeshData
    subgraph->inputs(wallBCOrchSM);

    // Three-phase pipeline: Orchestrator → Kernel (parallel) → Collector
    subgraph->edges(wallBCOrchSM, wallBCKernelTask);
    subgraph->edges(wallBCKernelTask, wallBCCollectorSM);

    // Exit: Collector emits MeshData
    subgraph->outputs(wallBCCollectorSM);

    return subgraph;
}

#endif // WALLBC_SUBGRAPH_H
