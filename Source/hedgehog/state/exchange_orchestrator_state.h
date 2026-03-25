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

/// Pure dependency gate for mesh exchange orchestration.
///
/// Receives MeshData "push done" signals from the upstream ExchangePushTask.
/// Each signal means that mesh NM has finished copying its data to all
/// same-rank targets.  The gate tracks which meshes have pushed and emits
/// a mesh downstream only when all of its receive-dependencies have also
/// pushed.
///
/// This state contains NO I/O — all copies and communication happen in
/// upstream/downstream tasks.  This makes it reusable across different
/// exchange codes (5, 3, 1, etc.) by pairing it with different push tasks.
///
/// A mesh can proceed when:
///   - It has been pushed (its data is saved in targets' OMESHes)
///   - All of its receive-dependencies have been pushed (its OMESH is fresh)
class ExchangeGateState
    : public hh::AbstractState<1, MeshData, MeshData> {
public:
    ExchangeGateState(std::shared_ptr<MeshDependencyGraph> depGraph)
        : depGraph_(std::move(depGraph)),
          satisfied_(static_cast<size_t>(depGraph_->totalMeshes())),
          lower_(depGraph_->lowerMesh()),
          upper_(depGraph_->upperMesh()) {

        int nmeshes = upper_ - lower_ + 1;
        pendingMeshes_.resize(static_cast<size_t>(depGraph_->totalMeshes() + 1));
        emitted_.resize(static_cast<size_t>(nmeshes), false);

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

        pendingMeshes_[static_cast<size_t>(nm)] = data;
        satisfied_.set(static_cast<size_t>(nm - 1));

        // Check all local meshes for newly-satisfied dependencies
        for (int destNM = lower_; destNM <= upper_; ++destNM) {
            size_t localIdx = static_cast<size_t>(destNM - lower_);
            if (emitted_[localIdx]) continue;
            if (!pendingMeshes_[static_cast<size_t>(destNM)]) continue;

            if (noDeps_[localIdx] ||
                satisfied_.containsAll(depGraph_->sameRankRecvDeps(destNM))) {
                emitted_[localIdx] = true;
                this->addResult(std::move(pendingMeshes_[static_cast<size_t>(destNM)]));
            }
        }

        auto t1 = std::chrono::steady_clock::now();
        gateTime_ += std::chrono::duration<double>(t1 - t0).count();
        ++invocations_;

        if (allEmitted()) {
            resetRound();
        }
    }

    [[nodiscard]] std::string info() const {
        std::ostringstream oss;
        oss << "EXCHANGE_GATE\\n"
            << std::fixed << std::setprecision(3)
            << "gate " << gateTime_ << "s"
            << " / " << invocations_ << " arrivals";
        if (invocations_ > 0) {
            oss << "\\navg " << std::setprecision(1)
                << (gateTime_ * 1e6 / invocations_) << "us";
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
    double gateTime_ = 0.0;
    int invocations_ = 0;
};

/// StateManager for ExchangeGateState, with dot-file diagnostics.
class ExchangeGateManager
    : public hh::StateManager<1, MeshData, MeshData> {
public:
    ExchangeGateManager(
        std::shared_ptr<ExchangeGateState> const &state,
        std::string const &name)
        : hh::StateManager<1, MeshData, MeshData>(
              state, name) {}

    [[nodiscard]] std::string extraPrintingInformation() const override {
        this->state()->lock();
        auto ret = std::dynamic_pointer_cast<ExchangeGateState>(
            this->state())->info();
        this->state()->unlock();
        return ret;
    }
};

#endif // EXCHANGE_ORCHESTRATOR_STATE_H
