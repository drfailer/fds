#ifndef MESH_DEPENDENCIES_MANAGER_STATE_H
#define MESH_DEPENDENCIES_MANAGER_STATE_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include <vector>
#include <chrono>
#include <iomanip>
#include <sstream>
#include "../data/mesh_data.h"
#include "../data/mesh_exchange_data.h"
#include "../tool/mesh_dependency_graph.h"

/// Dependency-aware mesh exchange orchestrator.
///
/// Receives MeshData tokens from the upstream kernel (one per mesh).
/// Tracks arrivals and emits MeshExchangeData tokens as soon as a mesh's
/// same-rank neighbors have all arrived (or been exchanged).
///
/// The emitted MeshExchangeData carries a list of neighbors that still
/// need bidirectional copies.  When mesh A is emitted and exchanged with
/// neighbor B, the state marks A as "done" in B's dependency tracking.
/// This can cascade: marking A done may satisfy B's dependencies, causing
/// B to be emitted in the same execute() call.
///
/// Single-process only (same-rank neighbors).  Cross-rank dependencies
/// will be handled by a future communicator task.
class MeshDependenciesManagerState
    : public hh::AbstractState<1, MeshData, MeshExchangeData> {
public:
    MeshDependenciesManagerState(std::shared_ptr<MeshDependencyGraph> depGraph)
        : depGraph_(std::move(depGraph)),
          lower_(depGraph_->lowerMesh()),
          upper_(depGraph_->upperMesh()) {

        int nmeshes = upper_ - lower_ + 1;
        int total = depGraph_->totalMeshes();

        // Per-mesh tracking (indexed by mesh number, 1-based)
        pendingMeshes_.resize(static_cast<size_t>(total + 1));
        arrived_.resize(static_cast<size_t>(nmeshes), false);
        emitted_.resize(static_cast<size_t>(nmeshes), false);

        // Per-mesh bitset: tracks which neighbors have arrived OR been exchanged
        // A mesh can be emitted when this bitset contains all its same-rank neighbors.
        satisfiedNeighbors_.resize(
            static_cast<size_t>(total + 1),
            DynBitset(static_cast<size_t>(total)));

        // Pre-compute which meshes have no same-rank neighbors (can be emitted immediately)
        noDeps_.resize(static_cast<size_t>(nmeshes), false);
        for (int nm = lower_; nm <= upper_; ++nm) {
            if (depGraph_->sameRankNeighbors(nm).count() == 0) {
                noDeps_[localIdx(nm)] = true;
            }
        }
    }

    void execute(std::shared_ptr<MeshData> data) override {
        auto t0 = std::chrono::steady_clock::now();

        int nm = data->nm;
        pendingMeshes_[static_cast<size_t>(nm)] = data;
        arrived_[localIdx(nm)] = true;

        // Mark nm as arrived in all local meshes' satisfied sets
        for (int destNM = lower_; destNM <= upper_; ++destNM) {
            if (depGraph_->sameRankNeighbors(destNM).contains(
                    static_cast<size_t>(nm - 1))) {
                satisfiedNeighbors_[static_cast<size_t>(destNM)].set(
                    static_cast<size_t>(nm - 1));
            }
        }

        // Cascade: check all local meshes for newly-satisfiable emissions
        emitReady();

        auto t1 = std::chrono::steady_clock::now();
        stateTime_ += std::chrono::duration<double>(t1 - t0).count();
        ++invocations_;

        if (allEmitted()) {
            resetRound();
        }
    }

    [[nodiscard]] std::string info() const {
        std::ostringstream oss;
        oss << "DEPS_MANAGER\\n"
            << std::fixed << std::setprecision(3)
            << "state " << stateTime_ << "s"
            << " / " << invocations_ << " arrivals";
        if (invocations_ > 0) {
            oss << "\\navg " << std::setprecision(1)
                << (stateTime_ * 1e6 / invocations_) << "us";
        }
        return oss.str();
    }

private:
    [[nodiscard]] size_t localIdx(int nm) const {
        return static_cast<size_t>(nm - lower_);
    }

    /// Try to emit all ready meshes, cascading as exchanges mark neighbors done.
    void emitReady() {
        bool progress = true;
        while (progress) {
            progress = false;
            for (int nm = lower_; nm <= upper_; ++nm) {
                size_t li = localIdx(nm);
                if (emitted_[li]) continue;
                if (!arrived_[li]) continue;

                // Check if all same-rank neighbors are satisfied
                if (!noDeps_[li] &&
                    !satisfiedNeighbors_[static_cast<size_t>(nm)].containsAll(
                        depGraph_->sameRankNeighbors(nm))) {
                    continue;
                }

                // Compute list of neighbors that haven't been exchanged yet
                std::vector<int> pendingNeighbors;
                for (int nom : depGraph_->sameRankNeighborsList(nm)) {
                    if (!emitted_[localIdx(nom)]) {
                        pendingNeighbors.push_back(nom);
                    }
                }

                emitted_[li] = true;

                // Mark nm as "exchanged" in all remaining local neighbors'
                // satisfied sets (this may trigger cascade emissions)
                for (int destNM = lower_; destNM <= upper_; ++destNM) {
                    if (depGraph_->sameRankNeighbors(destNM).contains(
                            static_cast<size_t>(nm - 1))) {
                        satisfiedNeighbors_[static_cast<size_t>(destNM)].set(
                            static_cast<size_t>(nm - 1));
                    }
                }

                this->addResult(std::make_shared<MeshExchangeData>(
                    std::move(pendingMeshes_[static_cast<size_t>(nm)]),
                    std::move(pendingNeighbors)));

                progress = true;  // re-scan for cascaded emissions
            }
        }
    }

    [[nodiscard]] bool allEmitted() const {
        for (size_t i = 0; i < emitted_.size(); ++i) {
            if (!emitted_[i]) return false;
        }
        return true;
    }

    void resetRound() {
        for (auto &p : pendingMeshes_) p.reset();
        std::fill(arrived_.begin(), arrived_.end(), false);
        std::fill(emitted_.begin(), emitted_.end(), false);
        for (int nm = lower_; nm <= upper_; ++nm) {
            satisfiedNeighbors_[static_cast<size_t>(nm)].reset();
        }
    }

    std::shared_ptr<MeshDependencyGraph> depGraph_;
    int lower_;
    int upper_;
    std::vector<std::shared_ptr<MeshData>> pendingMeshes_;
    std::vector<bool> arrived_;
    std::vector<bool> emitted_;
    std::vector<DynBitset> satisfiedNeighbors_;
    std::vector<bool> noDeps_;
    double stateTime_ = 0.0;
    int invocations_ = 0;
};

/// StateManager for MeshDependenciesManagerState, with dot-file diagnostics.
class MeshDependenciesManager
    : public hh::StateManager<1, MeshData, MeshExchangeData> {
public:
    MeshDependenciesManager(
        std::shared_ptr<MeshDependenciesManagerState> const &state,
        std::string const &name)
        : hh::StateManager<1, MeshData, MeshExchangeData>(
              state, name) {}

    [[nodiscard]] std::string extraPrintingInformation() const override {
        this->state()->lock();
        auto ret = std::dynamic_pointer_cast<MeshDependenciesManagerState>(
            this->state())->info();
        this->state()->unlock();
        return ret;
    }
};

#endif // MESH_DEPENDENCIES_MANAGER_STATE_H
