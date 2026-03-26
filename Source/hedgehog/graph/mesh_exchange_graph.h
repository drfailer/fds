#ifndef MESH_EXCHANGE_GRAPH_H
#define MESH_EXCHANGE_GRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_exchange_data.h"
#include "../data/termination_data.h"
#include "../state/mesh_dependencies_manager_state.h"
#include "../tool/mesh_dependency_graph.h"

/// Reusable sub-graph encapsulating the dependency-managed parallel exchange cycle.
///
/// Architecture:
///   MeshDepsManager(state) <-> ExchangeTask(task, parallel)
///
/// Input T tokens arrive, the dependency manager tracks per-mesh state
/// (NotArrived -> Wait -> Processing -> Processed -> Done) and dispatches
/// MeshExchangeData to the exchange task.  When all neighbors are processed,
/// T tokens are emitted to the output.
///
/// @tparam T Data type flowing through the exchange (must be accepted by
///           MeshDependenciesManagerState, typically MeshData)
template <typename T>
class MeshExchangeGraph : public hh::Graph<2, T, TerminationData, T> {
public:
    /// @param depGraph Pre-built mesh dependency graph
    /// @param exchangeTask Task performing the actual exchange copies
    MeshExchangeGraph(
        std::shared_ptr<MeshDependencyGraph> depGraph,
        std::shared_ptr<hh::AbstractTask<1, MeshExchangeData, MeshExchangeData>> exchangeTask)
        : hh::Graph<2, T, TerminationData, T>("MeshExchange") {

        auto depManagerSM = std::make_shared<MeshDependenciesManager>(
            std::make_shared<MeshDependenciesManagerState>(std::move(depGraph)),
            "MeshDepsManager");

        // Entry: T -> depManager (arrival tracking), TerminationData -> depManager (termination)
        this->template input<T>(depManagerSM);
        this->template input<TerminationData>(depManagerSM);

        // Cycle: depManager -> exchangeTask -> depManager
        this->edges(depManagerSM, exchangeTask);
        this->template edge<MeshExchangeData>(exchangeTask, depManagerSM);

        // Exit: depManager emits T (Done meshes)
        this->outputs(depManagerSM);
    }

    /// Wire TerminationData from a parent graph to this exchange sub-graph.
    ///
    /// Call this in the parent graph builder to route TerminationData:
    /// @code
    ///   MeshExchangeGraph<MeshData>::wireTermination(parentGraph, exchangeGraph);
    /// @endcode
    template <size_t S, typename... AllTypes>
    static void wireTermination(
        std::shared_ptr<hh::Graph<S, AllTypes...>> parentGraph,
        std::shared_ptr<MeshExchangeGraph<T>> self) {
        parentGraph->template input<TerminationData>(self);
    }
};

#endif // MESH_EXCHANGE_GRAPH_H
