#ifndef EXCHANGE_ORCHESTRATOR_STATE_H
#define EXCHANGE_ORCHESTRATOR_STATE_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include <vector>
#include <chrono>
#include <iomanip>
#include <sstream>
#include "../data/mesh_data.h"
#include "../tool/mesh_dependency_graph.h"
#include "../fds_fortran_interface.h"

/// Dependency-aware mesh exchange state (push-then-gate).
///
/// Replaces the global fds_mesh_exchange(5) barrier with per-mesh dependency
/// tracking.  When a mesh NM arrives:
///
///   1. PUSH: Copy NM's FVX/FVY/FVZ/H to all same-rank targets' OMESH(NM).
///      This reads NM's data BEFORE the pressure solve modifies it.
///
///   2. GATE: Mark NM as "pushed" in the satisfied_ bitset.  Check if any
///      local mesh's receive-dependencies are now fully satisfied.  For each
///      satisfied mesh, emit MeshData to the downstream pressure solve.
///
/// A mesh can proceed to the pressure solve only when:
///   - Its own data has been pushed (saved to targets' OMESHes)
///   - All of its receive-dependencies have pushed (its OMESH has fresh data)
///
/// Thread safety: The state is single-threaded (Hedgehog guarantee).  Push
/// copies run sequentially within execute(), which is correct because:
///   - Different pushes write to disjoint OMESH entries
///   - No mesh starts its solve before this state emits it
///   - The state reads source mesh data before any solve modifies it
class ExchangeOrchestratorState
    : public hh::AbstractState<1, MeshData, MeshData> {
public:
    ExchangeOrchestratorState(std::shared_ptr<MeshDependencyGraph> depGraph)
        : depGraph_(std::move(depGraph)),
          satisfied_(static_cast<size_t>(depGraph_->totalMeshes())),
          lower_(depGraph_->lowerMesh()),
          upper_(depGraph_->upperMesh()) {

        int nmeshes = upper_ - lower_ + 1;
        pendingMeshes_.resize(static_cast<size_t>(depGraph_->totalMeshes() + 1));
        emitted_.resize(static_cast<size_t>(nmeshes), false);

        // Pre-compute which local meshes have no same-rank dependencies
        noDeps_.resize(static_cast<size_t>(nmeshes), false);
        for (int nm = lower_; nm <= upper_; ++nm) {
            if (depGraph_->sameRankRecvDeps(nm).count() == 0) {
                noDeps_[static_cast<size_t>(nm - lower_)] = true;
            }
        }
    }

    void execute(std::shared_ptr<MeshData> data) override {
        auto t0 = std::chrono::steady_clock::now();

        int nm = data->nm;

        // Store the mesh pointer
        pendingMeshes_[static_cast<size_t>(nm)] = data;

        // PUSH: Copy nm's data to all same-rank targets' OMESHes.
        // Must happen NOW, before nm enters the pressure solve.
        for (int target : depGraph_->sendTargets(nm)) {
            if (fds_mesh_process(target) == fds_mesh_process(nm)) {
                fds_flux_copy_neighbor_ts(nm, target);
                ++copyCount_;
            }
        }

        // Mark nm as "pushed" (0-based bit index)
        satisfied_.set(static_cast<size_t>(nm - 1));

        // GATE: Check all local meshes for newly-satisfied dependencies
        for (int destNM = lower_; destNM <= upper_; ++destNM) {
            size_t localIdx = static_cast<size_t>(destNM - lower_);
            if (emitted_[localIdx]) continue;

            // destNM must have arrived (and pushed) itself
            if (!pendingMeshes_[static_cast<size_t>(destNM)]) continue;

            if (noDeps_[localIdx] ||
                satisfied_.containsAll(depGraph_->sameRankRecvDeps(destNM))) {
                emitted_[localIdx] = true;
                this->addResult(std::move(pendingMeshes_[static_cast<size_t>(destNM)]));
            }
        }

        auto t1 = std::chrono::steady_clock::now();
        orchTime_ += std::chrono::duration<double>(t1 - t0).count();
        ++invocations_;

        // Reset state when all local meshes have been emitted
        if (allEmitted()) {
            resetRound();
        }
    }

    [[nodiscard]] std::string info() const {
        std::ostringstream oss;
        oss << "EXCHANGE_ORCHESTRATOR\\n"
            << std::fixed << std::setprecision(3)
            << "total " << orchTime_ << "s"
            << " / " << invocations_ << " arrivals"
            << " / " << copyCount_ << " copies";
        if (invocations_ > 0) {
            oss << "\\navg " << std::setprecision(1)
                << (orchTime_ * 1e6 / invocations_) << "us/arrival";
        }
        return oss.str();
    }

private:
    [[nodiscard]] bool allEmitted() const {
        for (size_t i = 0; i < emitted_.size(); ++i) {
            if (!emitted_[i]) return false;
        }
        return true;
    }

    void resetRound() {
        satisfied_.reset();
        for (auto &p : pendingMeshes_) p.reset();
        std::fill(emitted_.begin(), emitted_.end(), false);
    }

    std::shared_ptr<MeshDependencyGraph> depGraph_;
    DynBitset satisfied_;
    int lower_;
    int upper_;
    std::vector<std::shared_ptr<MeshData>> pendingMeshes_;
    std::vector<bool> emitted_;
    std::vector<bool> noDeps_;
    double orchTime_ = 0.0;
    int invocations_ = 0;
    int copyCount_ = 0;
};

/// StateManager for ExchangeOrchestratorState, with dot-file diagnostics.
class ExchangeOrchestratorManager
    : public hh::StateManager<1, MeshData, MeshData> {
public:
    ExchangeOrchestratorManager(
        std::shared_ptr<ExchangeOrchestratorState> const &state,
        std::string const &name)
        : hh::StateManager<1, MeshData, MeshData>(
              state, name) {}

    [[nodiscard]] std::string extraPrintingInformation() const override {
        this->state()->lock();
        auto ret = std::dynamic_pointer_cast<ExchangeOrchestratorState>(
            this->state())->info();
        this->state()->unlock();
        return ret;
    }
};

#endif // EXCHANGE_ORCHESTRATOR_STATE_H
