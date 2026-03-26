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
#include "../data/termination_data.h"
#include "../tool/mesh_dependency_graph.h"

/// Per-mesh state within the exchange cycle.
enum class ExchangeState { NotArrived, Wait, Processing, Processed, Done };

/// Dependency-aware mesh exchange orchestrator with parallel pull-only exchange.
///
/// Manages a cycle with FluxExchangeTask:
///   - Receives MeshData from the upstream kernel (baroclinic)
///   - Dispatches MeshExchangeData to the exchange task when dependencies are met
///   - Receives MeshExchangeData back from the exchange task
///   - Emits MeshData to the downstream kernel (solve) when fully done
///
/// Per-mesh state machine:
///   NotArrived -> Wait (on arrival from upstream)
///   Wait -> Processing (when all neighbors have arrived)
///   Processing -> Processed (when exchange task returns)
///   Processed -> Done (when all neighbors are at least Processed)
///
/// Pull-only exchange: each mesh pulls data from ALL its same-rank neighbors.
/// Since each mesh writes only to its own OMESH buffers, different meshes
/// CAN be in Processing simultaneously without races.
///
/// The Done gate ensures no mesh proceeds to the solve task while any of its
/// neighbors is still in the exchange task.  This prevents a fast mesh from
/// starting solve (which modifies source data) while a slow neighbor is
/// still pulling from it.
///
/// Counter-based dependency tracking avoids full-mesh scans:
///   - unarrivedNeighborCount: guards Wait -> Processing
///   - unprocessedNeighborCount: guards Processed -> Done
///
/// Termination: receives TerminationData from the graph input when the
/// simulation is complete, setting done_=true so canTerminate() returns true.
class MeshDependenciesManagerState
    : public hh::AbstractState<3, MeshData, MeshExchangeData, TerminationData, MeshExchangeData, MeshData> {
public:
    MeshDependenciesManagerState(std::shared_ptr<MeshDependencyGraph> depGraph)
        : depGraph_(std::move(depGraph)),
          lower_(depGraph_->lowerMesh()),
          upper_(depGraph_->upperMesh()),
          nmeshes_(upper_ - lower_ + 1) {

        size_t n = static_cast<size_t>(nmeshes_);
        meshState_.resize(n, ExchangeState::NotArrived);
        pendingMeshes_.resize(n);
        unarrivedNeighborCount_.resize(n);
        unprocessedNeighborCount_.resize(n);

        for (int nm = lower_; nm <= upper_; ++nm) {
            size_t li = localIdx(nm);
            int nNeighbors = static_cast<int>(
                depGraph_->sameRankNeighborsList(nm).size());
            unarrivedNeighborCount_[li] = nNeighbors;
            unprocessedNeighborCount_[li] = nNeighbors;
        }
    }

    /// Handle arrival from upstream kernel (baroclinic).
    void execute(std::shared_ptr<MeshData> data) override {
        auto t0 = std::chrono::steady_clock::now();

        if (roundComplete_) {
            resetRound();
        }

        int nm = data->nm;
        size_t li = localIdx(nm);
        meshState_[li] = ExchangeState::Wait;
        pendingMeshes_[li] = data;

        for (int nom : depGraph_->sameRankNeighborsList(nm)) {
            unarrivedNeighborCount_[localIdx(nom)]--;
        }

        tryStartProcessing(nm);
        for (int nom : depGraph_->sameRankNeighborsList(nm)) {
            tryStartProcessing(nom);
        }

        auto t1 = std::chrono::steady_clock::now();
        stateTime_ += std::chrono::duration<double>(t1 - t0).count();
        ++arrivalCount_;
    }

    /// Handle return from exchange task (cycle back).
    void execute(std::shared_ptr<MeshExchangeData> data) override {
        auto t0 = std::chrono::steady_clock::now();

        int nm = data->mesh->nm;
        size_t li = localIdx(nm);
        meshState_[li] = ExchangeState::Processed;
        pendingMeshes_[li] = data->mesh;

        for (int nom : depGraph_->sameRankNeighborsList(nm)) {
            unprocessedNeighborCount_[localIdx(nom)]--;
        }

        tryFinish(nm);
        for (int nom : depGraph_->sameRankNeighborsList(nm)) {
            tryFinish(nom);
        }

        auto t1 = std::chrono::steady_clock::now();
        stateTime_ += std::chrono::duration<double>(t1 - t0).count();
        ++returnCount_;
    }

    /// Handle termination signal from graph input.
    void execute(std::shared_ptr<TerminationData>) override {
        done_ = true;
    }

    [[nodiscard]] bool isDone() const { return done_; }

    [[nodiscard]] std::string info() const {
        std::ostringstream oss;
        oss << "DEPS_MANAGER\\n"
            << std::fixed << std::setprecision(3)
            << "state " << stateTime_ << "s"
            << " / " << arrivalCount_ << " arrivals"
            << " / " << returnCount_ << " returns";
        int total = arrivalCount_ + returnCount_;
        if (total > 0) {
            oss << "\\navg " << std::setprecision(1)
                << (stateTime_ * 1e6 / total) << "us";
        }
        return oss.str();
    }

private:
    [[nodiscard]] size_t localIdx(int nm) const {
        return static_cast<size_t>(nm - lower_);
    }

    void resetRound() {
        std::fill(meshState_.begin(), meshState_.end(),
                  ExchangeState::NotArrived);
        for (auto &p : pendingMeshes_) p.reset();
        for (int nm = lower_; nm <= upper_; ++nm) {
            size_t li = localIdx(nm);
            int nNeighbors = static_cast<int>(
                depGraph_->sameRankNeighborsList(nm).size());
            unarrivedNeighborCount_[li] = nNeighbors;
            unprocessedNeighborCount_[li] = nNeighbors;
        }
        doneCount_ = 0;
        roundComplete_ = false;
    }

    /// Attempt Wait -> Processing transition.
    /// Condition: all same-rank neighbors have arrived.
    void tryStartProcessing(int nm) {
        size_t li = localIdx(nm);
        if (meshState_[li] != ExchangeState::Wait) return;
        if (unarrivedNeighborCount_[li] > 0) return;

        meshState_[li] = ExchangeState::Processing;

        // Pull-only: each mesh pulls from ALL its same-rank neighbors.
        this->addResult(std::make_shared<MeshExchangeData>(
            std::move(pendingMeshes_[li]),
            depGraph_->sameRankNeighborsList(nm)));
    }

    /// Attempt Processed -> Done transition.
    /// Condition: all same-rank neighbors are at least Processed.
    void tryFinish(int nm) {
        size_t li = localIdx(nm);
        if (meshState_[li] != ExchangeState::Processed) return;
        if (unprocessedNeighborCount_[li] > 0) return;

        meshState_[li] = ExchangeState::Done;
        ++doneCount_;

        this->addResult(std::move(pendingMeshes_[li]));

        if (doneCount_ == nmeshes_) {
            roundComplete_ = true;
        }
    }

    std::shared_ptr<MeshDependencyGraph> depGraph_;
    int lower_;
    int upper_;
    int nmeshes_;

    std::vector<ExchangeState> meshState_;
    std::vector<std::shared_ptr<MeshData>> pendingMeshes_;
    std::vector<int> unarrivedNeighborCount_;
    std::vector<int> unprocessedNeighborCount_;

    int doneCount_ = 0;
    bool roundComplete_ = false;
    bool done_ = false;

    double stateTime_ = 0.0;
    int arrivalCount_ = 0;
    int returnCount_ = 0;
};

/// StateManager for MeshDependenciesManagerState with cycle termination support.
class MeshDependenciesManager
    : public hh::StateManager<3, MeshData, MeshExchangeData, TerminationData, MeshExchangeData, MeshData> {
public:
    MeshDependenciesManager(
        std::shared_ptr<MeshDependenciesManagerState> const &state,
        std::string const &name)
        : hh::StateManager<3, MeshData, MeshExchangeData, TerminationData, MeshExchangeData, MeshData>(
              state, name) {}

    [[nodiscard]] bool canTerminate() const override {
        this->state()->lock();
        auto s = std::dynamic_pointer_cast<MeshDependenciesManagerState>(
            this->state());
        bool ret = s->isDone();
        this->state()->unlock();
        return ret;
    }

    [[nodiscard]] std::string extraPrintingInformation() const override {
        this->state()->lock();
        auto ret = std::dynamic_pointer_cast<MeshDependenciesManagerState>(
            this->state())->info();
        this->state()->unlock();
        return ret;
    }
};

#endif // MESH_DEPENDENCIES_MANAGER_STATE_H
