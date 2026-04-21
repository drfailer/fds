#ifndef COPY_TASK_H
#define COPY_TASK_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "exchange_bundle.h"
#include "mesh_dependency_graph.h"
#include "exchange_dispatch.h"
#include "../fds_fortran_interface.h"

/// Variadic copy task: performs intra-rank copies and cross-rank unpacks
/// for all exchange types in a single node.
template<MeshState... Ks>
class CopyTask : public BundleDispatch<
    CopyTask<Ks...>,
    hh::AbstractTask<sizeof...(Ks), ExchangeBundle<Ks>..., MeshData<Ks>...>,
    Ks...>
{
    using HHBase = hh::AbstractTask<sizeof...(Ks), ExchangeBundle<Ks>..., MeshData<Ks>...>;
    using DispatchBase = BundleDispatch<CopyTask<Ks...>, HHBase, Ks...>;
    friend DispatchBase;

public:
    CopyTask(size_t numThreads, std::shared_ptr<MeshDependencyGraph> depGraph)
        : DispatchBase("Copy", numThreads),
          depGraph_(std::move(depGraph)) {}

    template<MeshState K>
    void handleBundle(std::shared_ptr<ExchangeBundle<K>> bundle) {
        int nm = bundle->mesh->nm;

        constexpr int code = PackData<K>::exchangeCode();
        for (int nom : depGraph_->sameRankNeighborsList(nm)) {
            fds_exchange_copy_neighbor_ts(code, nom, nm);
        }

        for (auto &pd : bundle->remoteHalos) {
            pd->unpackIntoOMesh();
        }

        this->addResult(std::move(bundle->mesh));
    }

    std::shared_ptr<HHBase> copy() override {
        return std::make_shared<CopyTask<Ks...>>(this->numberThreads(), depGraph_);
    }

private:
    std::shared_ptr<MeshDependencyGraph> depGraph_;
};

#endif // COPY_TASK_H
