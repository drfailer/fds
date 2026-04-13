#ifndef MESH_EXCHANGE_GRAPH_H
#define MESH_EXCHANGE_GRAPH_H

#include <hedgehog/hedgehog.h>
#include <hedgehog_comm.h>
#include <memory>
#include <algorithm>
#include <array>
#include <vector>
#include <unordered_map>
#include "../data/mesh_data.h"
#include "../data/exchange_mesh_data.h"
#include "../tool/mesh_dependency_graph.h"
#include "../tool/exchange_buffer.h"
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
/// FanOutTask pushes slab data into ExchangeMeshData objects.
/// CommunicatorTask routes: same-rank loop-back, cross-rank via MPI.
/// WriteBufferTask copies slab data into ExchangeBuffer, emits DepSignal.
/// Gate waits for pending token + all dep signals, then emits downstream.
/// PullTask reads from ExchangeBuffer into OMESH arrays.
///
/// Double-buffered: 2 ExchangeBuffer copies per code, indexed by roundId % 2.
/// The Gate is also double-buffered (2 RoundSlots). This allows concurrent
/// exchange rounds (e.g. pre-solve round 0 and post-solve round 1) to
/// proceed without data races.
///
/// The exchange code is selected at runtime via md->exchangeCode.
class MeshExchangeGraph : public hh::Graph<1, MeshData, MeshData> {
    using Pool = hh::comm::tool::MemoryPool<ExchangeMeshData>;
    using BufferMap = std::unordered_map<int, std::array<std::shared_ptr<ExchangeBuffer>, 2>>;
public:
    /// @param depGraph Pre-built mesh dependency graph
    /// @param pushThreads Thread count for the fan-out task
    /// @param pullThreads Thread count for the write-buffer and pull tasks
    /// @param supportedCodes Exchange codes this graph handles (e.g. {5})
    /// @param commService Optional MPI comm service; when non-null a
    ///                    CommunicatorTask is wired between fan-out and write
    /// @param name Graph name for profiling/debugging
    MeshExchangeGraph(
        std::shared_ptr<MeshDependencyGraph> depGraph,
        size_t pushThreads,
        size_t pullThreads,
        std::vector<int> const &supportedCodes,
        hh::comm::CommService *commService = nullptr,
        std::string const &name = "MeshExchange")
        : hh::Graph<1, MeshData, MeshData>(name) {

        // Build double-buffered per-code ExchangeBuffers
        auto buffers = std::make_shared<BufferMap>();
        int globalMax = 0;
        for (int code : supportedCodes) {
            auto buf0 = std::make_shared<ExchangeBuffer>(*depGraph, code);
            auto buf1 = std::make_shared<ExchangeBuffer>(*depGraph, code);
            buffers->emplace(code, std::array<std::shared_ptr<ExchangeBuffer>, 2>{buf0, buf1});
            int codeMax = fds_exchange_max_slab_size(code);
            globalMax = std::max(globalMax, codeMax);
        }
        ExchangeMeshData::globalMaxSlabSize_ = globalMax;

        // Shared pool for ExchangeMeshData lifecycle
        auto pool = std::make_shared<Pool>();
        int totalTargets = computeTotalSendTargets(*depGraph);
        pool->template fill<ExchangeMeshData>(
            static_cast<size_t>(std::max(4 * totalTargets, 32)));

        auto fanOutTask = std::make_shared<ExchangeFanOutTask>(
            pushThreads, depGraph, pool);
        auto writeTask = std::make_shared<ExchangeWriteBufferTask>(
            pullThreads, buffers, pool);
        auto gateSM = std::make_shared<ExchangeDepsGateManager>(
            std::make_shared<ExchangeDepsGateState>(depGraph),
            name + "_Gate");
        auto pullTask = std::make_shared<ExchangePullBufferTask>(
            pullThreads, buffers);

        // Graph input broadcasts MeshData to BOTH fanOutTask and gate
        this->template input<MeshData>(fanOutTask);
        this->template input<MeshData>(gateSM);

        if (commService) {
            auto commTask = std::make_shared<hh::CommunicatorTask<ExchangeMeshData>>(
                commService, name + "_Comm");
            commTask->template strategy<ExchangeMeshData>(
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
