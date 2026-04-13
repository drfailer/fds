#ifndef EXCHANGE_DEPS_STATE_H
#define EXCHANGE_DEPS_STATE_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include <vector>
#include <array>
#include <chrono>
#include <iomanip>
#include <sstream>
#include "../data/mesh_data.h"
#include "../data/exchange_mesh_data.h"
#include "../tool/mesh_dependency_graph.h"

/// Double-buffered dependency gate state for the mesh exchange pipeline.
///
/// Receives two input types:
///   - MeshData: pending tokens from graph input (one per local mesh per round)
///   - ExchangeDepSignal: dep satisfaction from WriteBufferTask
///     (one per (source, dest) pair per round)
///
/// Uses roundId/exchangeRound % 2 to index into one of two RoundSlots,
/// so concurrent exchange rounds (e.g. pre-solve round 0 and post-solve
/// round 1) never share state.
///
/// When all recv deps for a local mesh are satisfied in a slot (or the mesh
/// has no deps), emits its pending MeshData token downstream to the pull task.
///
/// Slot resets automatically after all local meshes have been emitted from
/// that slot.
class ExchangeDepsGateState
    : public hh::AbstractState<2, MeshData, ExchangeDepSignal, MeshData> {
public:
    ExchangeDepsGateState(std::shared_ptr<MeshDependencyGraph> depGraph)
        : depGraph_(std::move(depGraph)),
          lower_(depGraph_->lowerMesh()),
          upper_(depGraph_->upperMesh()),
          nmeshes_(upper_ - lower_ + 1) {

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

        // Initialize both round slots
        size_t totalPlus1 = static_cast<size_t>(depGraph_->totalMeshes() + 1);
        for (auto &slot : slots_) {
            slot.pendingMeshes.resize(totalPlus1);
            slot.arrivalCount.resize(static_cast<size_t>(nmeshes_), 0);
            slot.emitCount = 0;
        }
    }

    /// Pending token from graph input broadcast.
    void execute(std::shared_ptr<MeshData> data) override {
        auto t0 = std::chrono::steady_clock::now();
        int slot = data->exchangeRound % 2;
        slots_[slot].pendingMeshes[static_cast<size_t>(data->nm)] = data;
        tryEmit(data->nm, slot);
        auto t1 = std::chrono::steady_clock::now();
        gateTime_ += std::chrono::duration<double>(t1 - t0).count();
        ++invocations_;
    }

    /// Dep arrival: sourceNm wrote data for destNom into buffer.
    void execute(std::shared_ptr<ExchangeDepSignal> signal) override {
        auto t0 = std::chrono::steady_clock::now();
        int slot = signal->roundId % 2;
        int idx = signal->destNom - lower_;
        if (idx >= 0 && idx < nmeshes_) {
            slots_[slot].arrivalCount[static_cast<size_t>(idx)]++;
            tryEmit(signal->destNom, slot);
        }
        auto t1 = std::chrono::steady_clock::now();
        gateTime_ += std::chrono::duration<double>(t1 - t0).count();
        ++invocations_;
    }

    [[nodiscard]] std::string info() const {
        std::ostringstream oss;
        oss << "EXCHANGE_DEPS_GATE (2-slot)\\n"
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
    struct RoundSlot {
        std::vector<std::shared_ptr<MeshData>> pendingMeshes;
        std::vector<int> arrivalCount;
        int emitCount = 0;
    };

    void tryEmit(int nm, int slot) {
        int idx = nm - lower_;
        if (idx < 0 || idx >= nmeshes_) return;
        size_t uidx = static_cast<size_t>(idx);
        auto &s = slots_[slot];
        if (!s.pendingMeshes[static_cast<size_t>(nm)]) return;
        if (noDeps_[uidx] || s.arrivalCount[uidx] == totalRecvDeps_[uidx]) {
            this->addResult(
                std::move(s.pendingMeshes[static_cast<size_t>(nm)]));
            if (++s.emitCount == nmeshes_) {
                resetSlot(slot);
            }
        }
    }

    void resetSlot(int slot) {
        auto &s = slots_[slot];
        for (auto &p : s.pendingMeshes) p.reset();
        std::fill(s.arrivalCount.begin(), s.arrivalCount.end(), 0);
        s.emitCount = 0;
    }

    std::shared_ptr<MeshDependencyGraph> depGraph_;
    int lower_;
    int upper_;
    int nmeshes_;
    std::array<RoundSlot, 2> slots_;
    std::vector<int> totalRecvDeps_;  // read-only, shared across slots
    std::vector<bool> noDeps_;        // read-only, shared across slots
    double gateTime_ = 0.0;
    int invocations_ = 0;
};

/// StateManager wrapper for ExchangeDepsGateState.
/// No canTerminate needed — the exchange graph has no internal cycle.
/// Cycle-breaking is handled by the PostExchangeRouter in the parent graph.
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
