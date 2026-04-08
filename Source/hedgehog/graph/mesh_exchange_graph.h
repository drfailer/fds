#ifndef MESH_EXCHANGE_GRAPH_H
#define MESH_EXCHANGE_GRAPH_H

#include <hedgehog/hedgehog.h>
#include <communicator_task.hpp>
#include <memory>
#include "../data/mesh_data.h"
#include "../tool/mesh_dependency_graph.h"
#include "../tool/exchange_buffer.h"
#include "../tool/exchange_strategy.h"
#include "../task/exchange_push_buffer_task.h"
#include "../task/exchange_pull_buffer_task.h"
#include "../state/exchange_deps_state.h"

/// Reusable sub-graph encapsulating the push/buffer/pull exchange pipeline.
///
/// Architecture (linear — no internal cycle):
///   [CommunicatorTask →] PushTask(parallel) → DepsGate(1 thread) → PullTask(parallel)
///
/// When a CommService is provided, a CommunicatorTask sits at the front of
/// the pipeline.  For same-rank data (strategy returns {self}), the
/// communicator calls addResult directly (no MPI).  For cross-rank data
/// (future Phase 3), the communicator serializes and sends via MPI.
///
/// Push copies source mesh arrays into a pre-allocated ExchangeBuffer.
/// DepsGate emits downstream only when all receive dependencies have pushed.
/// Pull copies from the ExchangeBuffer into OMESH arrays.
///
/// Each instance owns its own ExchangeBuffer, so concurrent exchanges
/// (e.g. pre-solve and post-solve in the pressure iteration) never share
/// buffers and cannot race.
///
/// No TerminationData needed — there is no internal cycle.  The parent
/// graph's cycle management (e.g. pressure convergence barrier) handles
/// termination.
///
/// @tparam T Data type flowing through the exchange (typically MeshData)
/// @tparam Strategy Exchange strategy trait (e.g. FluxExchangeStrategy)
template <typename T, typename Strategy = FluxExchangeStrategy>
class MeshExchangeGraph : public hh::Graph<1, T, T> {
public:
    /// @param depGraph Pre-built mesh dependency graph
    /// @param pushThreads Thread count for the push task
    /// @param pullThreads Thread count for the pull task
    /// @param commService Optional MPI comm service; when non-null a CommunicatorTask
    ///                    is wired at the front of the pipeline
    /// @param name Graph name for profiling/debugging
    MeshExchangeGraph(
        std::shared_ptr<MeshDependencyGraph> depGraph,
        size_t pushThreads,
        size_t pullThreads,
        hh::comm::CommService *commService = nullptr,
        std::string const &name = "MeshExchange")
        : hh::Graph<1, T, T>(name) {

        auto buffer = std::make_shared<ExchangeBuffer<Strategy>>(*depGraph);

        auto pushTask = std::make_shared<ExchangePushBufferTask<Strategy>>(
            pushThreads, depGraph, buffer);
        auto gateTask = std::make_shared<ExchangeDepsGateTask>(
            std::move(depGraph));
        auto pullTask = std::make_shared<ExchangePullBufferTask<Strategy>>(
            pullThreads, buffer);

        if (commService) {
            // CommunicatorTask at the front: loop-back strategy (same rank)
            auto commTask = std::make_shared<hh::CommunicatorTask<T>>(
                commService, name + "_Comm");
            commTask->template strategy<T>(
                [rank = commService->rank()](auto) {
                    return std::vector<hh::comm::rank_t>{rank};
                });

            this->template input<T>(commTask);
            this->edges(commTask, pushTask);
        } else {
            this->template input<T>(pushTask);
        }

        this->edges(pushTask, gateTask);
        this->edges(gateTask, pullTask);
        this->outputs(pullTask);
    }
};

#endif // MESH_EXCHANGE_GRAPH_H
