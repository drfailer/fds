#ifndef EXCHANGE_DEPS_STATE_H
#define EXCHANGE_DEPS_STATE_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include <vector>
#include <chrono>
#include <iomanip>
#include <sstream>
#include "../data/mesh_data.h"
#include "../data/exchange_flux_mesh_data.h"
#include "../tool/mesh_dependency_graph.h"

/// Dependency gate state for the mesh exchange pipeline.
///
/// Receives two input types:
///   - MeshData: pending tokens from graph input (one per local mesh per round)
///   - ExchangeDepSignal: dep satisfaction from WriteBufferTask
///     (one per (source, dest) pair)
///
/// When all recv deps for a local mesh are satisfied (or the mesh has no deps),
/// emits its pending MeshData token downstream to the pull task.
///
/// Round resets automatically after all local meshes have been emitted.
class ExchangeDepsGateState
    : public hh::AbstractState<2, MeshData, ExchangeDepSignal, MeshData> {
public:
    ExchangeDepsGateState(std::shared_ptr<MeshDependencyGraph> depGraph)
        : depGraph_(std::move(depGraph)),
          lower_(depGraph_->lowerMesh()),
          upper_(depGraph_->upperMesh()),
          nmeshes_(upper_ - lower_ + 1) {

        pendingMeshes_.resize(static_cast<size_t>(depGraph_->totalMeshes() + 1));
        arrivalCount_.resize(static_cast<size_t>(nmeshes_), 0);
        totalRecvDeps_.resize(static_cast<size_t>(nmeshes_), 0);
        noDeps_.resize(static_cast<size_t>(nmeshes_), false);

        for (int nm = lower_; nm <= upper_; ++nm) {
            size_t idx = static_cast<size_t>(nm - lower_);
            int depCount = static_cast<int>(depGraph_->recvDeps(nm).count());
            totalRecvDeps_[idx] = depCount;
            if (depCount == 0) {
                noDeps_[idx] = true;
            }
        }
    }

    /// Pending token from graph input broadcast.
    void execute(std::shared_ptr<MeshData> data) override {
        auto t0 = std::chrono::steady_clock::now();
        pendingMeshes_[static_cast<size_t>(data->nm)] = data;
        tryEmit(data->nm);
        auto t1 = std::chrono::steady_clock::now();
        gateTime_ += std::chrono::duration<double>(t1 - t0).count();
        ++invocations_;
    }

    /// Dep arrival: sourceNm wrote data for destNom into buffer.
    void execute(std::shared_ptr<ExchangeDepSignal> signal) override {
        auto t0 = std::chrono::steady_clock::now();
        int idx = signal->destNom - lower_;
        if (idx >= 0 && idx < nmeshes_) {
            arrivalCount_[static_cast<size_t>(idx)]++;
            tryEmit(signal->destNom);
        }
        auto t1 = std::chrono::steady_clock::now();
        gateTime_ += std::chrono::duration<double>(t1 - t0).count();
        ++invocations_;
    }

    [[nodiscard]] std::string info() const {
        std::ostringstream oss;
        oss << "EXCHANGE_DEPS_GATE\\n"
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
    void tryEmit(int nm) {
        int idx = nm - lower_;
        if (idx < 0 || idx >= nmeshes_) return;
        size_t uidx = static_cast<size_t>(idx);
        if (!pendingMeshes_[static_cast<size_t>(nm)]) return;
        if (noDeps_[uidx] || arrivalCount_[uidx] == totalRecvDeps_[uidx]) {
            this->addResult(
                std::move(pendingMeshes_[static_cast<size_t>(nm)]));
            emitCount_++;
            if (emitCount_ == nmeshes_) {
                resetRound();
            }
        }
    }

    void resetRound() {
        for (auto &p : pendingMeshes_) p.reset();
        std::fill(arrivalCount_.begin(), arrivalCount_.end(), 0);
        emitCount_ = 0;
    }

    std::shared_ptr<MeshDependencyGraph> depGraph_;
    int lower_;
    int upper_;
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> pendingMeshes_;
    std::vector<int> arrivalCount_;
    std::vector<int> totalRecvDeps_;
    std::vector<bool> noDeps_;
    int emitCount_ = 0;
    double gateTime_ = 0.0;
    int invocations_ = 0;
};

/// StateManager wrapper for ExchangeDepsGateState.
/// No canTerminate needed — there is no internal cycle in the exchange graph.
class ExchangeDepsGateManager
    : public hh::StateManager<2, MeshData, ExchangeDepSignal, MeshData> {
public:
    ExchangeDepsGateManager(
        std::shared_ptr<ExchangeDepsGateState> const &state,
        std::string const &name)
        : hh::StateManager<2, MeshData, ExchangeDepSignal, MeshData>(
              state, name) {}

    [[nodiscard]] std::string extraPrintingInformation() const override {
        this->state()->lock();
        auto ret = std::dynamic_pointer_cast<ExchangeDepsGateState>(
            this->state())->info();
        this->state()->unlock();
        return ret;
    }
};

#endif // EXCHANGE_DEPS_STATE_H
