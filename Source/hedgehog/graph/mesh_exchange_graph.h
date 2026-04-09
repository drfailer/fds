#ifndef MESH_EXCHANGE_GRAPH_H
#define MESH_EXCHANGE_GRAPH_H

#include <hedgehog/hedgehog.h>
#include <hedgehog_comm.h>
#include <memory>
#include <algorithm>
#include "../data/mesh_data.h"
#include "../data/exchange_flux_mesh_data.h"
#include "../tool/mesh_dependency_graph.h"
#include "../tool/exchange_buffer.h"
#include "../tool/exchange_strategy.h"
#include "../task/exchange_fanout_task.h"
#include "../task/exchange_write_buffer_task.h"
#include "../task/exchange_pull_buffer_task.h"
#include "../state/exchange_deps_state.h"

/// Reusable sub-graph encapsulating the mesh exchange pipeline.
///
/// Architecture:
///   Input -> FanOutTask(parallel) -> [CommTask] -> WriteBufferTask(parallel)
///                                                        |DepSignal
///   Input -> Gate(StateManager) <------------------------+
///              |
///              +-> PullTask(parallel) -> Output
///
/// Graph input broadcasts MeshData to both FanOutTask and Gate.
/// FanOutTask pushes slab data into ExchangeFluxMeshData objects.
/// CommunicatorTask routes: same-rank loop-back, cross-rank via MPI.
/// WriteBufferTask copies slab data into ExchangeBuffer, emits DepSignal.
/// Gate waits for pending token + all dep signals, then emits downstream.
/// PullTask reads from ExchangeBuffer into OMESH arrays.
///
/// Each instance owns its own ExchangeBuffer, so concurrent exchanges
/// (e.g. pre-solve and post-solve in the pressure iteration) never share
/// buffers and cannot race.
///
/// @tparam T Data type flowing through the exchange (typically MeshData)
/// @tparam Strategy Exchange strategy trait (e.g. FluxExchangeStrategy)
template <typename T, typename Strategy = FluxExchangeStrategy>
class MeshExchangeGraph : public hh::Graph<1, T, T> {
    using EFD = ExchangeFluxMeshData<Strategy>;
    using Pool = hh::comm::tool::MemoryPool<EFD>;
public:
    /// @param depGraph Pre-built mesh dependency graph
    /// @param pushThreads Thread count for the fan-out task
    /// @param pullThreads Thread count for the write-buffer and pull tasks
    /// @param commService Optional MPI comm service; when non-null a
    ///                    CommunicatorTask is wired between fan-out and write
    /// @param name Graph name for profiling/debugging
    MeshExchangeGraph(
        std::shared_ptr<MeshDependencyGraph> depGraph,
        size_t pushThreads,
        size_t pullThreads,
        hh::comm::CommService *commService = nullptr,
        std::string const &name = "MeshExchange")
        : hh::Graph<1, T, T>(name) {

        auto buffer = std::make_shared<ExchangeBuffer<Strategy>>(*depGraph);

        // Shared pool for ExchangeFluxMeshData lifecycle
        auto pool = std::make_shared<Pool>();
        int totalTargets = computeTotalSendTargets(*depGraph);
        pool->template fill<EFD>(
            static_cast<size_t>(std::max(4 * totalTargets, 32)));

        auto fanOutTask = std::make_shared<ExchangeFanOutTask<Strategy>>(
            pushThreads, depGraph, pool);
        auto writeTask = std::make_shared<ExchangeWriteBufferTask<Strategy>>(
            pullThreads, buffer, pool);
        auto gateSM = std::make_shared<ExchangeDepsGateManager>(
            std::make_shared<ExchangeDepsGateState>(depGraph),
            name + "_Gate");
        auto pullTask = std::make_shared<ExchangePullBufferTask<Strategy>>(
            pullThreads, buffer);

        // Graph input broadcasts MeshData to BOTH fanOutTask and gate
        this->template input<T>(fanOutTask);
        this->template input<T>(gateSM);

        if (commService) {
            auto commTask = std::make_shared<hh::CommunicatorTask<EFD>>(
                commService, name + "_Comm");
            commTask->template strategy<EFD>(
                [](auto data) {
                    return std::vector<hh::comm::rank_t>{
                        static_cast<hh::comm::rank_t>(
                            fds_mesh_process(data->header_.destNom))};
                });
            commTask->setMemoryManager(pool);

            this->edges(fanOutTask, commTask);
            this->edges(commTask, writeTask);
        } else {
            this->edges(fanOutTask, writeTask);
        }

        this->template edge<ExchangeDepSignal>(writeTask, gateSM);
        this->edges(gateSM, pullTask);
        this->outputs(pullTask);
    }

private:
    static int computeTotalSendTargets(const MeshDependencyGraph &depGraph) {
        int total = 0;
        for (int nm = depGraph.lowerMesh(); nm <= depGraph.upperMesh(); ++nm) {
            total += static_cast<int>(depGraph.sendTargets(nm).size());
        }
        return total;
    }
};

#endif // MESH_EXCHANGE_GRAPH_H
