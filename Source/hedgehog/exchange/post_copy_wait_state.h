#ifndef POST_COPY_WAIT_STATE_H
#define POST_COPY_WAIT_STATE_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include <vector>
#include <sstream>
#include "../data/mesh_data.h"
#include "../data/termination_data.h"
#include "mesh_dependency_graph.h"
#include "pack_data.h"
#include "exchange_dispatch.h"

/// Variadic post-copy gate: ensures no mesh proceeds downstream while a
/// neighbor's CopyTask is still reading from its arrays.
///
/// Same-rank-only tracking. On release, retags MeshData<K> → MeshData<Out>
/// using the compile-time K→Out mapping from ExchKind pairs.
///
/// Single-threaded sequential task (1 thread). Receives TerminationData for
/// cycle termination — canTerminate returns isDone directly, no Manager needed.
template<typename... ExchTypes>
class PostCopyWaitTask;

template<MeshState... Ks, MeshState... Outs>
class PostCopyWaitTask<ExchKind<Ks, Outs>...>
    : public MeshDataDispatch<
          PostCopyWaitTask<ExchKind<Ks, Outs>...>,
          hh::AbstractTask<sizeof...(Ks) + 1,
              MeshData<Ks>..., TerminationData,
              MeshData<Outs>...>,
          Ks...>
{
    using HHBase = hh::AbstractTask<sizeof...(Ks) + 1,
        MeshData<Ks>..., TerminationData,
        MeshData<Outs>...>;
    using DispatchBase = MeshDataDispatch<
        PostCopyWaitTask<ExchKind<Ks, Outs>...>, HHBase, Ks...>;
    friend DispatchBase;

    template<MeshState K>
    struct PerTypeData {
        std::vector<bool> arrived;
        std::vector<std::shared_ptr<MeshData<K>>> pendingMeshes;
        std::vector<int> unarrivedCount;
        int doneCount = 0;
        bool roundComplete = false;
    };

    template<MeshState K>
    static constexpr MeshState outputTagFor() {
        constexpr MeshState ks[] = {Ks...};
        constexpr MeshState outs[] = {Outs...};
        for (size_t i = 0; i < sizeof...(Ks); ++i) {
            if (ks[i] == K) return outs[i];
        }
        return K;
    }

public:
    explicit PostCopyWaitTask(std::shared_ptr<MeshDependencyGraph> depGraph)
        : DispatchBase("PostCopyWait", 1),
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
        pt.arrived[li] = true;
        pt.pendingMeshes[li] = data;

        for (int nom : depGraph_->sameRankNeighborsList(nm)) {
            pt.unarrivedCount[localIdx(nom)]--;
        }

        tryRelease<K>(nm);
        for (int nom : depGraph_->sameRankNeighborsList(nm)) {
            tryRelease<K>(nom);
        }
    }

private:
    [[nodiscard]] size_t localIdx(int nm) const {
        return static_cast<size_t>(nm - lower_);
    }

    template<MeshState K>
    void initType() {
        auto &pt = std::get<PerTypeData<K>>(perType_);
        size_t n = static_cast<size_t>(nmeshes_);
        pt.arrived.resize(n, false);
        pt.pendingMeshes.resize(n);
        pt.unarrivedCount.resize(n);
    }

    template<MeshState K>
    void resetType() {
        auto &pt = std::get<PerTypeData<K>>(perType_);
        std::fill(pt.arrived.begin(), pt.arrived.end(), false);
        for (auto &p : pt.pendingMeshes) p.reset();
        for (int nm = lower_; nm <= upper_; ++nm) {
            pt.unarrivedCount[localIdx(nm)] = static_cast<int>(
                depGraph_->sameRankNeighborsList(nm).size());
        }
        pt.doneCount = 0;
        pt.roundComplete = false;
    }

    template<MeshState K>
    void tryRelease(int nm) {
        auto &pt = std::get<PerTypeData<K>>(perType_);
        size_t li = localIdx(nm);
        if (!pt.arrived[li]) return;
        if (pt.unarrivedCount[li] > 0) return;

        pt.arrived[li] = false;
        ++pt.doneCount;

        constexpr MeshState Out = outputTagFor<K>();
        this->addResult(retag<Out>(std::move(pt.pendingMeshes[li])));

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

#endif // POST_COPY_WAIT_STATE_H
