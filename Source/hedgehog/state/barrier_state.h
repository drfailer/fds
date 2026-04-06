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
#include "../data/barrier_data.h"
#include "../fds_fortran_interface.h"

/// Barrier function signature: takes collected meshes, runs global work.
using BarrierFn = std::function<void(std::vector<std::shared_ptr<MeshData>>&)>;

/// Barrier collector task that collects N MeshData tokens, runs a barrier
/// function, then re-emits N MeshData tokens.
///
/// Runs on a single thread. Accumulates tokens internally, firing the
/// barrier function only when all expected tokens have arrived.
///
/// Timing and routine info are exposed via extraPrintingInformation()
/// for dot-file diagnostics.
class BarrierCollectorTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    /// @param nmeshes Number of unique meshes (determines output count)
    /// @param name Task name for dot-file display
    /// @param routines Label for dot-file display
    /// @param fn Barrier function called when all tokens arrive
    /// @param totalExpected Total tokens to collect before firing (default: nmeshes).
    ///        Set to numBranches*nmeshes for fork-join+barrier merges.
    BarrierCollectorTask(int nmeshes, std::string name, std::string routines,
                         BarrierFn fn, int totalExpected = 0)
        : hh::AbstractTask<1, MeshData, MeshData>(std::move(name), 1),
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
            this->batchAddResult(collected_);
            for (auto &md : collected_) { md = nullptr; }
        }
    }

    [[nodiscard]] std::string extraPrintingInformation() const override {
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

/// Helper to create a BarrierCollectorTask.
inline auto makeBarrierSM(int nmeshes, std::string name,
                           std::string routines, BarrierFn fn,
                           int totalExpected = 0) {
    return std::make_shared<BarrierCollectorTask>(
        nmeshes, std::move(name), std::move(routines), std::move(fn),
        totalExpected);
}

/// Task variant: accepts pre-collected BarrierData, runs barrier function,
/// scatters N MeshData tokens downstream. Used when upstream already holds
/// collected meshes (e.g. retry loop output).
class BarrierTask : public hh::AbstractTask<1, BarrierData, MeshData> {
public:
    BarrierTask(std::string name, std::string routines, BarrierFn fn)
        : hh::AbstractTask<1, BarrierData, MeshData>(std::move(name), 1),
          routines_(std::move(routines)), fn_(std::move(fn)) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        auto t0 = std::chrono::steady_clock::now();
        fn_(data->meshes);
        auto t1 = std::chrono::steady_clock::now();
        totalTime_ += std::chrono::duration<double>(t1 - t0).count();
        ++invocations_;
        this->batchAddResult(data->meshes);
    }

    [[nodiscard]] std::string extraPrintingInformation() const override {
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
    std::string routines_;
    BarrierFn fn_;
    double totalTime_ = 0.0;
    int invocations_ = 0;
};

/// Helper to create a BarrierTask (BarrierData → MeshData).
inline auto makeBarrierTask(std::string name, std::string routines, BarrierFn fn) {
    return std::make_shared<BarrierTask>(std::move(name), std::move(routines), std::move(fn));
}

#endif // BARRIER_STATE_H
