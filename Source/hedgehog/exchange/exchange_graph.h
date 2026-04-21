#ifndef EXCHANGE_GRAPH_H
#define EXCHANGE_GRAPH_H

#include <hedgehog/hedgehog.h>
#include <hedgehog_comm.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/termination_data.h"
#include "mesh_dependency_graph.h"
#include "exchange_dispatch.h"
#include "pack_task.h"
#include "pack_data.h"
#include "copy_task.h"
#include "neighbor_wait_state.h"
#include "post_copy_wait_state.h"

/// Unified dependency-aware exchange graph handling multiple exchange types
/// through a single set of nodes.
///
/// Template: ExchangeGraph<ExchKind<K1,Out1>, ExchKind<K2,Out2>, ...>
///
/// One PackTask, one NeighborWait, one CopyTask, one PostCopyWait — each
/// handles all exchange types via variadic templates and tuple storage.
/// One CommunicatorTask handles all MPI exchanges.
///
/// Multi-process (commService != nullptr):
///
///   MeshData<Ki> -> PackTask -> PackData<Ki> -> Comm --+
///                      |                               |
///                  MeshData<Ki>                    PackData<Ki>
///                      |                               |
///                      +---> NeighborWait <-------------+
///                                 |
///                          ExchangeBundle<Ki>
///                                 |
///                             CopyTask
///                                 |
///                          MeshData<Ki>
///                                 |
///                          PostCopyWait
///                                 |
///                          MeshData<Outi> -> output
///
/// Single-process (commService == nullptr):
///   MeshData<Ki> -> NeighborWait -> CopyTask -> PostCopyWait -> output
template<typename... ExchTypes>
class ExchangeGraph;

template<MeshState... Ks, MeshState... Outs>
class ExchangeGraph<ExchKind<Ks, Outs>...> : public hh::Graph<
    sizeof...(Ks) + 1,
    MeshData<Ks>..., TerminationData,
    MeshData<Outs>...>
{
    using GraphBase = hh::Graph<
        sizeof...(Ks) + 1,
        MeshData<Ks>..., TerminationData,
        MeshData<Outs>...>;

public:
    ExchangeGraph(
        std::shared_ptr<MeshDependencyGraph> depGraph,
        hh::comm::CommService *commService = nullptr,
        std::string const &name = "Exchange")
        : GraphBase(name) {

        auto waitTask = std::make_shared<NeighborWaitTask<Ks...>>(depGraph);

        auto copyTask = std::make_shared<CopyTask<Ks...>>(1, depGraph);

        auto postCopyTask = std::make_shared<PostCopyWaitTask<ExchKind<Ks, Outs>...>>(depGraph);

        if (commService) {
            auto packTask = std::make_shared<PackTask<Ks...>>(depGraph);

            auto commTask = std::make_shared<hh::CommunicatorTask<PackData<Ks>...>>(
                commService, name + "_Comm");

            (commTask->template strategy<PackData<Ks>>([](auto data) {
                return std::vector<hh::comm::rank_t>{
                    static_cast<hh::comm::rank_t>(
                        fds_mesh_process(data->destNom))};
            }), ...);

            size_t poolSize = static_cast<size_t>(depGraph->totalMeshes()) * 2;
            memoryPool_ = std::make_shared<hh::comm::tool::MemoryPool<PackData<Ks>...>>();
            (memoryPool_->template fill<PackData<Ks>>(poolSize), ...);
            commTask->setMemoryManager(memoryPool_);

            // Input -> PackTask
            (this->template input<MeshData<Ks>>(packTask), ...);

            // PackTask -> NeighborWait (MeshData pass-through)
            (this->template edge<MeshData<Ks>>(packTask, waitTask), ...);

            // PackTask -> Comm (PackData to send)
            (this->template edge<PackData<Ks>>(packTask, commTask), ...);

            // Comm -> NeighborWait (PackData received from remote)
            (this->template edge<PackData<Ks>>(commTask, waitTask), ...);
        } else {
            (this->template input<MeshData<Ks>>(waitTask), ...);
        }

        // NeighborWait -> CopyTask (ExchangeBundle)
        (this->template edge<ExchangeBundle<Ks>>(waitTask, copyTask), ...);

        // CopyTask -> PostCopyWait (MeshData)
        (this->template edge<MeshData<Ks>>(copyTask, postCopyTask), ...);

        // TerminationData -> both wait tasks
        this->template input<TerminationData>(waitTask);
        this->template input<TerminationData>(postCopyTask);

        this->outputs(postCopyTask);
    }

private:
    std::shared_ptr<hh::comm::tool::MemoryPool<PackData<Ks>...>> memoryPool_;
};

#endif // EXCHANGE_GRAPH_H
