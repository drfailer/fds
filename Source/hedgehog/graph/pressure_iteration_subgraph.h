#ifndef PRESSURE_ITERATION_SUBGRAPH_H
#define PRESSURE_ITERATION_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include <communicator_task.hpp>
#include <tool/memory_pool.hpp>
#include "../data/pressure_iteration_data.h"
#include "../data/mesh_data.h"
#include "../data/flux_exchange_data.h"
#include "../data/termination_signal.h"
#include "../task/pressure_iteration_tasks.h"
#include "../state/pressure_iteration_state.h"
#include "../state/flux_exchange_state.h"

/// Build the pressure iteration sub-graph.
///
/// Architecture with per-neighbor flux exchange:
///
///   PreCollector(state) -> FluxPack(state) -> FluxCollector(state)
///                                               -> SolveKernel(task)
///                                                    -> PostLoop(state)
///                ^                                          | (cycle)
///                +------------------------------------------+
///
/// Single-process: FluxPack emits FluxExchangeData and MeshData directly to
/// FluxCollector (no CommunicatorTask needed — all data is same-rank).
///
/// Multi-process (future): a CommunicatorTask is inserted between FluxPack
/// and FluxCollector to route FluxExchangeData across MPI ranks.
///
/// FluxCollectorState emits MeshData when a mesh has received flux data
/// from all its neighbors (per-neighbor barrier, not global barrier).
///
/// @param tEnd Simulation end time
/// @param nmeshes Number of local meshes
/// @param kernelThreads Number of threads for parallel kernel tasks
/// @param predictor True for predictor phase
/// @param termSignal Shared termination signal
/// @param commService Pointer to the MPI comm service for the communicator task
/// @param presFlag Pressure solver flag (FFT_FLAG=0, ULMAT_FLAG=3)
inline auto buildPressureIterationSubgraph(double tEnd, int nmeshes,
                                            size_t kernelThreads,
                                            bool predictor,
                                            std::shared_ptr<TerminationSignal> termSignal,
                                            hh::comm::CommService *commService,
                                            int presFlag = 0) {
    using SubGraphType = hh::Graph<1, MeshData, MeshData>;
    auto subgraph = std::make_shared<SubGraphType>("PressureIteration");

    int nmOffset = fds_get_lower_mesh_index();

    // --- Phase 1: PreCollector (init + baroclinic, no exchange) ---
    auto preCollSM = std::make_shared<PressurePreCollectorManager>(
        std::make_shared<PressurePreCollector>(nmeshes),
        "PressurePreCollector");

    // --- Flux exchange: Pack -> Collect ---
    auto fluxPackSM = std::make_shared<FluxPackStateManager>(
        std::make_shared<FluxPackState>(nmeshes, nmOffset),
        "FluxPack");

    auto fluxCollectSM = std::make_shared<FluxCollectorStateManager>(
        std::make_shared<FluxCollectorState>(nmeshes, nmOffset),
        "FluxCollect");

    // --- Phase 2: Parallel pressure solve kernel ---
    auto solveKernel = std::make_shared<PressureSolveKernelTask>(kernelThreads, presFlag);

    // --- Phase 3: PostLoop (exchange(5) + velocity error + convergence) ---
    auto postLoopSM = std::make_shared<PressurePostLoopManager>(
        std::make_shared<PressurePostLoop>(nmeshes, tEnd, predictor, termSignal),
        "PressurePostLoop");

    // --- Wire the sub-graph ---

    // Entry: MeshData -> PreCollector
    subgraph->inputs(preCollSM);

    // PreCollector -> FluxPack (MeshData)
    subgraph->edges(preCollSM, fluxPackSM);

    bool multiProcess = commService && commService->nbProcesses() > 1;

    if (multiProcess) {
        // Multi-process: route FluxExchangeData through CommunicatorTask
        auto commTask = std::make_shared<hh::CommunicatorTask<FluxExchangeData>>(
            commService, "FluxExchange");
        // Send strategy: compute destination rank from destNM
        // For now: same rank (TODO: compute PROCESS(destNM) for cross-rank)
        commTask->strategy<FluxExchangeData>([commService](auto data) {
            return std::vector<hh::comm::rank_t>{commService->rank()};
        });
        auto mm = std::make_shared<hh::comm::tool::MemoryPool<FluxExchangeData>>();
        mm->fill<FluxExchangeData>(2);
        commTask->setMemoryManager(mm);

        subgraph->edge<FluxExchangeData>(fluxPackSM, commTask);
        subgraph->edge<MeshData>(fluxPackSM, fluxCollectSM);
        subgraph->edge<FluxExchangeData>(commTask, fluxCollectSM);
    } else {
        // Single-process: FluxPack -> FluxCollect directly (no MPI needed)
        subgraph->edges(fluxPackSM, fluxCollectSM);
    }

    // FluxCollector -> SolveKernel (MeshData)
    subgraph->edges(fluxCollectSM, solveKernel);

    // SolveKernel -> PostLoop (MeshData)
    subgraph->edges(solveKernel, postLoopSM);

    // Cycle: PressureIterData -> back to PreCollector
    subgraph->edge<PressureIterData>(postLoopSM, preCollSM);

    // Exit: MeshData -> subgraph output
    subgraph->outputs(postLoopSM);

    return subgraph;
}

#endif // PRESSURE_ITERATION_SUBGRAPH_H
