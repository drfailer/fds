#ifndef BARRIER_STATE_H
#define BARRIER_STATE_H

#include <hedgehog/hedgehog.h>
#include <chrono>
#include <functional>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Generic barrier state that collects N MeshData tokens, runs a barrier
/// function, then re-emits N MeshData tokens.  Replaces the common pattern
/// of CollectorState(MeshData->BarrierData) + BarrierTask(BarrierData->MeshData)
/// with a single state node, eliminating one inter-node queue.
///
/// Timing and routine info are exposed via the BarrierStateManager's
/// extraPrintingInformation() for dot-file diagnostics.
using BarrierFn = std::function<void(std::vector<std::shared_ptr<MeshData>>&)>;

class BarrierState : public hh::AbstractState<1, MeshData, MeshData> {
public:
    /// @param nmeshes Number of unique meshes (determines output count)
    /// @param routines Label for dot-file display
    /// @param fn Barrier function called when all tokens arrive
    /// @param totalExpected Total tokens to collect before firing (default: nmeshes).
    ///        Set to numBranches*nmeshes for fork-join+barrier merges.
    BarrierState(int nmeshes, std::string routines, BarrierFn fn,
                 int totalExpected = 0)
        : hh::AbstractState<1, MeshData, MeshData>(),
          nmeshes_(nmeshes), routines_(std::move(routines)), fn_(std::move(fn)),
          totalExpected_(totalExpected > 0 ? totalExpected : nmeshes),
          nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_[data->nm - nmOffset_] = data;
        if (++count_ == totalExpected_) {
            auto t0 = std::chrono::steady_clock::now();
            fn_(collected_);
            auto t1 = std::chrono::steady_clock::now();
            totalTime_ += std::chrono::duration<double>(t1 - t0).count();
            ++invocations_;
            count_ = 0;
            for (auto &md : collected_) {
                this->addResult(md);
                md = nullptr;
            }
        }
    }

    [[nodiscard]] std::string info() const {
        std::ostringstream oss;
        oss << routines_ << "\\n"
            << std::fixed << std::setprecision(3) << totalTime_ << "s"
            << " / " << invocations_ << " calls";
        if (invocations_ > 0)
            oss << " / avg " << std::setprecision(3)
                << (totalTime_ * 1000.0 / invocations_) << "ms";
        return oss.str();
    }

private:
    int nmeshes_, nmOffset_, count_ = 0, totalExpected_;
    std::string routines_;
    BarrierFn fn_;
    double totalTime_ = 0.0;
    int invocations_ = 0;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

/// StateManager wrapping BarrierState, with extraPrintingInformation()
/// delegating to the state's info() method.
class BarrierStateManager
    : public hh::StateManager<1, MeshData, MeshData> {
public:
    BarrierStateManager(std::shared_ptr<BarrierState> const &state,
                        std::string const &name)
        : hh::StateManager<1, MeshData, MeshData>(state, name) {}

    [[nodiscard]] std::string extraPrintingInformation() const override {
        this->state()->lock();
        auto ret = std::dynamic_pointer_cast<BarrierState>(
            this->state())->info();
        this->state()->unlock();
        return ret;
    }
};

/// Helper to create a BarrierStateManager wrapping a BarrierState.
inline auto makeBarrierSM(int nmeshes, std::string name,
                           std::string routines, BarrierFn fn,
                           int totalExpected = 0) {
    return std::make_shared<BarrierStateManager>(
        std::make_shared<BarrierState>(nmeshes, std::move(routines), std::move(fn),
                                       totalExpected),
        std::move(name));
}

#endif // BARRIER_STATE_H
