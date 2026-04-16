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
/// Templated on MeshState so barriers can work with any MeshData variant.
template<MeshState S = MeshState::Default>
using BarrierFn = std::function<void(std::vector<std::shared_ptr<MeshData<S>>>&)>;

/// Barrier collector task that collects N MeshData tokens, runs a barrier
/// function, then re-emits N MeshData tokens.
///
/// Runs on a single thread. Accumulates tokens internally, firing the
/// barrier function only when all expected tokens have arrived.
template<MeshState S = MeshState::Default>
class BarrierCollectorTask
    : public hh::AbstractTask<1, MeshData<S>, MeshData<S>> {
public:
    BarrierCollectorTask(int nmeshes, std::string name, std::string routines,
                         BarrierFn<S> fn, int totalExpected = 0)
        : hh::AbstractTask<1, MeshData<S>, MeshData<S>>(std::move(name), 1),
          nmeshes_(nmeshes), routines_(std::move(routines)), fn_(std::move(fn)),
          totalExpected_(totalExpected > 0 ? totalExpected : nmeshes),
          nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData<S>> data) override {
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

    std::shared_ptr<hh::AbstractTask<1, MeshData<S>, MeshData<S>>> copy() override {
        return std::make_shared<BarrierCollectorTask<S>>(
            nmeshes_, routines_, routines_, fn_, totalExpected_);
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
    BarrierFn<S> fn_;
    double totalTime_ = 0.0;
    int invocations_ = 0;
    std::vector<std::shared_ptr<MeshData<S>>> collected_;
};

/// Helper to create a BarrierCollectorTask.
/// Uses generic Fn parameter so callers can pass auto-lambdas without
/// breaking template argument deduction on S.
template<MeshState S = MeshState::Default, typename Fn>
inline auto makeBarrierSM(int nmeshes, std::string name,
                           std::string routines, Fn&& fn,
                           int totalExpected = 0) {
    return std::make_shared<BarrierCollectorTask<S>>(
        nmeshes, std::move(name), std::move(routines),
        std::forward<Fn>(fn), totalExpected);
}

/// Task variant: accepts pre-collected BarrierData, runs barrier function,
/// scatters N MeshData tokens downstream. Used when upstream already holds
/// collected meshes (e.g. retry loop output).
class BarrierTask : public hh::AbstractTask<1, BarrierData, MeshData<>> {
public:
    BarrierTask(std::string name, std::string routines, BarrierFn<> fn)
        : hh::AbstractTask<1, BarrierData, MeshData<>>(std::move(name), 1),
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
    BarrierFn<> fn_;
    double totalTime_ = 0.0;
    int invocations_ = 0;
};

/// Helper to create a BarrierTask (BarrierData -> MeshData).
inline auto makeBarrierTask(std::string name, std::string routines, BarrierFn<> fn) {
    return std::make_shared<BarrierTask>(std::move(name), std::move(routines), std::move(fn));
}

/// Chain task: accepts BarrierData, runs barrier function, passes BarrierData
/// through WITHOUT scattering to MeshData.  Use between sequential barrier
/// operations to avoid unnecessary scatter-then-collect overhead.
class BarrierChainTask : public hh::AbstractTask<1, BarrierData, BarrierData> {
public:
    BarrierChainTask(std::string name, std::string routines, BarrierFn<> fn)
        : hh::AbstractTask<1, BarrierData, BarrierData>(std::move(name), 1),
          routines_(std::move(routines)), fn_(std::move(fn)) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        auto t0 = std::chrono::steady_clock::now();
        fn_(data->meshes);
        auto t1 = std::chrono::steady_clock::now();
        totalTime_ += std::chrono::duration<double>(t1 - t0).count();
        ++invocations_;
        this->addResult(data);
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
    BarrierFn<> fn_;
    double totalTime_ = 0.0;
    int invocations_ = 0;
};

/// Helper to create a BarrierChainTask (BarrierData -> BarrierData).
inline auto makeBarrierChainTask(std::string name, std::string routines, BarrierFn<> fn) {
    return std::make_shared<BarrierChainTask>(std::move(name), std::move(routines), std::move(fn));
}

/// Dual-input barrier: collects N MeshData<> + 1 BarrierData, runs barrier
/// function on the N meshes, then scatters N MeshData<> downstream.
/// The BarrierData input is consumed for counting only (meshes discarded).
/// Used when one fork branch produces BarrierData (e.g. radiation collector)
/// and the other produces N MeshData tokens.
class DualInputBarrierTask
    : public hh::AbstractTask<2, MeshData<>, BarrierData, MeshData<>> {
public:
    template<typename Fn>
    DualInputBarrierTask(int nmeshes, std::string name, std::string routines, Fn&& fn)
        : hh::AbstractTask<2, MeshData<>, BarrierData, MeshData<>>(
              std::move(name), 1),
          nmeshes_(nmeshes), routines_(std::move(routines)),
          fn_(std::forward<Fn>(fn)),
          nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData<>> data) override {
        collected_[data->nm - nmOffset_] = data;
        if (++count_ == nmeshes_ + 1) fire();
    }

    void execute(std::shared_ptr<BarrierData>) override {
        if (++count_ == nmeshes_ + 1) fire();
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
    void fire() {
        auto t0 = std::chrono::steady_clock::now();
        fn_(collected_);
        auto t1 = std::chrono::steady_clock::now();
        totalTime_ += std::chrono::duration<double>(t1 - t0).count();
        ++invocations_;
        count_ = 0;
        this->batchAddResult(collected_);
        for (auto &md : collected_) { md = nullptr; }
    }

    int nmeshes_, nmOffset_, count_ = 0;
    std::string routines_;
    BarrierFn<> fn_;
    double totalTime_ = 0.0;
    int invocations_ = 0;
    std::vector<std::shared_ptr<MeshData<>>> collected_;
};

/// Helper to create a DualInputBarrierTask.
template<typename Fn>
inline auto makeDualInputBarrier(int nmeshes, std::string name,
                                  std::string routines, Fn&& fn) {
    return std::make_shared<DualInputBarrierTask>(
        nmeshes, std::move(name), std::move(routines), std::forward<Fn>(fn));
}

/// Barrier collector that collects N MeshData, runs a barrier function,
/// and emits a single BarrierData (not N MeshData). Use when the barrier
/// performs global work and downstream only needs a completion signal.
/// The meshes are stored in the emitted BarrierData for downstream access.
class BarrierCollectToOneTask
    : public hh::AbstractTask<1, MeshData<>, BarrierData> {
public:
    template<typename Fn>
    BarrierCollectToOneTask(int nmeshes, std::string name, std::string routines, Fn&& fn)
        : hh::AbstractTask<1, MeshData<>, BarrierData>(std::move(name), 1),
          nmeshes_(nmeshes), routines_(std::move(routines)),
          fn_(std::forward<Fn>(fn)),
          nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData<>> data) override {
        collected_[data->nm - nmOffset_] = data;
        if (++count_ == nmeshes_) {
            auto t0 = std::chrono::steady_clock::now();
            fn_(collected_);
            auto t1 = std::chrono::steady_clock::now();
            totalTime_ += std::chrono::duration<double>(t1 - t0).count();
            ++invocations_;
            count_ = 0;
            auto bd = std::make_shared<BarrierData>();
            bd->meshes = std::move(collected_);
            collected_.resize(nmeshes_, nullptr);
            this->addResult(bd);
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
    int nmeshes_, nmOffset_, count_ = 0;
    std::string routines_;
    BarrierFn<> fn_;
    double totalTime_ = 0.0;
    int invocations_ = 0;
    std::vector<std::shared_ptr<MeshData<>>> collected_;
};

/// Helper to create a BarrierCollectToOneTask.
template<typename Fn>
inline auto makeBarrierCollectToOne(int nmeshes, std::string name,
                                     std::string routines, Fn&& fn) {
    return std::make_shared<BarrierCollectToOneTask>(
        nmeshes, std::move(name), std::move(routines), std::forward<Fn>(fn));
}

#endif // BARRIER_STATE_H
