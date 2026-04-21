#ifndef PACK_TASK_H
#define PACK_TASK_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "pack_data.h"
#include "mesh_dependency_graph.h"
#include "exchange_dispatch.h"

/// Variadic pack task: packs halo data for all exchange types in a single node.
///
/// For each MeshData<Ki> received, iterates remote neighbors and emits
/// PackData<Ki> (to CommunicatorTask) + MeshData<Ki> pass-through (to wait).
template<MeshState... Ks>
class PackTask : public MeshDataDispatch<
    PackTask<Ks...>,
    hh::AbstractTask<sizeof...(Ks), MeshData<Ks>..., MeshData<Ks>..., PackData<Ks>...>,
    Ks...>
{
    using HHBase = hh::AbstractTask<sizeof...(Ks), MeshData<Ks>..., MeshData<Ks>..., PackData<Ks>...>;
    using DispatchBase = MeshDataDispatch<PackTask<Ks...>, HHBase, Ks...>;
    friend DispatchBase;

public:
    PackTask(std::shared_ptr<MeshDependencyGraph> depGraph)
        : DispatchBase("Pack", 1),
          depGraph_(std::move(depGraph)) {}

    template<MeshState K>
    void handleMeshData(std::shared_ptr<MeshData<K>> md) {
        int nm = md->nm;
        for (int nom : depGraph_->sendTargets(nm)) {
            if (fds_mesh_process(nom) != fds_mesh_process(nm)) {
                auto pd = std::make_shared<PackData<K>>();
                pd->reset(nm, nom);
                pd->fillBuffer();
                this->addResult(pd);
            }
        }
        this->addResult(md);
    }

    std::shared_ptr<HHBase> copy() override {
        return std::make_shared<PackTask<Ks...>>(depGraph_);
    }

private:
    std::shared_ptr<MeshDependencyGraph> depGraph_;
};

#endif // PACK_TASK_H
