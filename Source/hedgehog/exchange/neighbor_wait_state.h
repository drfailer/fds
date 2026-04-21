#ifndef NEIGHBOR_WAIT_STATE_H
#define NEIGHBOR_WAIT_STATE_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include <vector>
#include <sstream>
#include "../data/mesh_data.h"
#include "../data/termination_data.h"
#include "mesh_dependency_graph.h"
#include "pack_data.h"
#include "exchange_bundle.h"
#include "exchange_dispatch.h"

/// Variadic dependency-aware gate tracking mesh and halo arrivals for all
/// exchange types in a single node.
///
/// Per-type tracking is stored in a tuple of PerTypeData<K> structs, accessed
/// via std::get<PerTypeData<K>>(perType_). Each type's rounds are independent.
///
/// Single-threaded sequential task (1 thread). Receives TerminationData for
/// cycle termination — canTerminate returns isDone directly, no Manager needed.
template<MeshState... Ks>
class NeighborWaitTask;

template<MeshState... Ks>
using NWHHBase_ = hh::AbstractTask<
    sizeof...(Ks) * 2 + 1,
    MeshData<Ks>..., PackData<Ks>..., TerminationData,
    ExchangeBundle<Ks>...>;

template<MeshState... Ks>
using NWMeshLayer_ = MeshDataDispatch<NeighborWaitTask<Ks...>, NWHHBase_<Ks...>, Ks...>;

template<MeshState... Ks>
using NWDispatchBase_ = PackDataDispatch<NeighborWaitTask<Ks...>, NWMeshLayer_<Ks...>, Ks...>;

template<MeshState... Ks>
class NeighborWaitTask : public NWDispatchBase_<Ks...> {
    template<typename, typename, MeshState...> friend struct MeshDataDispatch;
    template<typename, typename, MeshState...> friend struct PackDataDispatch;

    template<MeshState K>
    struct PerTypeData {
        std::vector<bool> localArrived;
        std::vector<std::shared_ptr<MeshData<K>>> pendingMeshes;
        std::vector<std::vector<std::shared_ptr<PackData<K>>>> pendingRemote;
        std::vector<int> unarrivedCount;
        int doneCount = 0;
        bool roundComplete = false;
    };

public:
    explicit NeighborWaitTask(std::shared_ptr<MeshDependencyGraph> depGraph)
        : NWDispatchBase_<Ks...>("NeighborWait", 1),
          depGraph_(std::move(depGraph)),
          lower_(depGraph_->lowerMesh()),
          upper_(depGraph_->upperMesh()),
          nmeshes_(upper_ - lower_ + 1) {
        (initType<Ks>(), ...);
        (resetType<Ks>(), ...);
    }

    void execute(std::shared_ptr<TerminationData>) override { done_ = true; }

    [[nodiscard]] bool canTerminate() const override { return done_; }

    template<MeshState K>
    void handleMeshData(std::shared_ptr<MeshData<K>> data) {
        auto &pt = std::get<PerTypeData<K>>(perType_);
        if (pt.roundComplete) resetType<K>();

        int nm = data->nm;
        size_t li = localIdx(nm);
        pt.localArrived[li] = true;
        pt.pendingMeshes[li] = data;

        for (int nom : depGraph_->sameRankNeighborsList(nm)) {
            pt.unarrivedCount[localIdx(nom)]--;
        }

        tryRelease<K>(nm);
        for (int nom : depGraph_->sameRankNeighborsList(nm)) {
            tryRelease<K>(nom);
        }
    }

    template<MeshState K>
    void handlePackData(std::shared_ptr<PackData<K>> pd) {
        auto &pt = std::get<PerTypeData<K>>(perType_);
        if (pt.roundComplete) resetType<K>();

        int destNom = pd->destNom;
        size_t li = localIdx(destNom);

        pt.pendingRemote[li].push_back(pd);
        pt.unarrivedCount[li]--;
        tryRelease<K>(destNom);
    }

private:
    [[nodiscard]] size_t localIdx(int nm) const {
        return static_cast<size_t>(nm - lower_);
    }

    template<MeshState K>
    void initType() {
        auto &pt = std::get<PerTypeData<K>>(perType_);
        size_t n = static_cast<size_t>(nmeshes_);
        pt.localArrived.resize(n, false);
        pt.pendingMeshes.resize(n);
        pt.pendingRemote.resize(n);
        pt.unarrivedCount.resize(n);
    }

    template<MeshState K>
    void resetType() {
        auto &pt = std::get<PerTypeData<K>>(perType_);
        std::fill(pt.localArrived.begin(), pt.localArrived.end(), false);
        for (auto &p : pt.pendingMeshes) p.reset();
        for (auto &v : pt.pendingRemote) v.clear();
        for (int nm = lower_; nm <= upper_; ++nm) {
            size_t li = localIdx(nm);
            int sameRank = static_cast<int>(
                depGraph_->sameRankNeighborsList(nm).size());
            int totalRecv = static_cast<int>(
                depGraph_->recvDeps(nm).count());
            int crossRank = totalRecv - static_cast<int>(
                depGraph_->sameRankRecvDeps(nm).count());
            pt.unarrivedCount[li] = sameRank + crossRank;
        }
        pt.doneCount = 0;
        pt.roundComplete = false;
    }

    template<MeshState K>
    void tryRelease(int nm) {
        auto &pt = std::get<PerTypeData<K>>(perType_);
        size_t li = localIdx(nm);
        if (!pt.localArrived[li]) return;
        if (pt.unarrivedCount[li] > 0) return;

        pt.localArrived[li] = false;
        ++pt.doneCount;

        auto bundle = std::make_shared<ExchangeBundle<K>>(
            std::move(pt.pendingMeshes[li]));
        bundle->remoteHalos = std::move(pt.pendingRemote[li]);
        this->addResult(std::move(bundle));

        if (pt.doneCount == nmeshes_) {
            pt.roundComplete = true;
        }
    }

    std::shared_ptr<MeshDependencyGraph> depGraph_;
    int lower_;
    int upper_;
    int nmeshes_;
    std::tuple<PerTypeData<Ks>...> perType_;
    bool done_ = false;
};

#endif // NEIGHBOR_WAIT_STATE_H
