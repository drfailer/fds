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
#include "../data/termination_data.h"
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

/// Eager dual-input barrier: collects N MeshData<> + 1 BarrierData.
/// Runs barrier function as soon as N MeshData arrive (doesn't wait for BarrierData).
/// Emits N MeshData only when BOTH the barrier function has run AND BarrierData arrived.
class EagerDualInputBarrierTask
    : public hh::AbstractTask<2, MeshData<>, BarrierData, MeshData<>> {
public:
    template<typename Fn>
    EagerDualInputBarrierTask(int nmeshes, std::string name, std::string routines, Fn&& fn)
        : hh::AbstractTask<2, MeshData<>, BarrierData, MeshData<>>(
              std::move(name), 1),
          nmeshes_(nmeshes), routines_(std::move(routines)),
          fn_(std::forward<Fn>(fn)),
          nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData<>> data) override {
        collected_[data->nm - nmOffset_] = data;
        if (++meshCount_ == nmeshes_) {
            auto t0 = std::chrono::steady_clock::now();
            fn_(collected_);
            auto t1 = std::chrono::steady_clock::now();
            totalTime_ += std::chrono::duration<double>(t1 - t0).count();
            ++invocations_;
            fnDone_ = true;
            tryEmit();
        }
    }

    void execute(std::shared_ptr<BarrierData>) override {
        exchDone_ = true;
        tryEmit();
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
    void tryEmit() {
        if (fnDone_ && exchDone_) {
            fnDone_ = false;
            exchDone_ = false;
            meshCount_ = 0;
            this->batchAddResult(collected_);
            for (auto &md : collected_) { md = nullptr; }
        }
    }

    int nmeshes_, nmOffset_, meshCount_ = 0;
    bool fnDone_ = false, exchDone_ = false;
    std::string routines_;
    BarrierFn<> fn_;
    double totalTime_ = 0.0;
    int invocations_ = 0;
    std::vector<std::shared_ptr<MeshData<>>> collected_;
};

/// Helper to create an EagerDualInputBarrierTask.
template<typename Fn>
inline auto makeEagerDualInputBarrier(int nmeshes, std::string name,
                                       std::string routines, Fn&& fn) {
    return std::make_shared<EagerDualInputBarrierTask>(
        nmeshes, std::move(name), std::move(routines), std::forward<Fn>(fn));
}

/// Barrier collector that collects N MeshData<S>, runs a barrier function,
/// and emits a single BarrierData (not N MeshData). Use when the barrier
/// performs global work and downstream only needs a completion signal.
/// The meshes are stored in the emitted BarrierData for downstream access.
/// Templated on MeshState so upstream tasks can use typed routing.
/// Internally retags to Default for the barrier function and BarrierData storage.
template<MeshState S = MeshState::Default>
class BarrierCollectToOneTask
    : public hh::AbstractTask<1, MeshData<S>, BarrierData> {
public:
    template<typename Fn>
    BarrierCollectToOneTask(int nmeshes, std::string name, std::string routines, Fn&& fn)
        : hh::AbstractTask<1, MeshData<S>, BarrierData>(std::move(name), 1),
          nmeshes_(nmeshes), routines_(std::move(routines)),
          fn_(std::forward<Fn>(fn)),
          nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData<S>> data) override {
        collected_[data->nm - nmOffset_] = retag<MeshState::Default>(data);
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
template<MeshState S = MeshState::Default, typename Fn>
inline auto makeBarrierCollectToOne(int nmeshes, std::string name,
                                     std::string routines, Fn&& fn) {
    return std::make_shared<BarrierCollectToOneTask<S>>(
        nmeshes, std::move(name), std::move(routines), std::forward<Fn>(fn));
}

/// Retagging barrier: collects N MeshData<InS>, runs barrier function, emits
/// N MeshData<OutS>. Used when adjacent pipeline stages use different MeshState
/// tags for type-based routing (e.g. packed parallel task ↔ sequential barriers).
template<MeshState InS, MeshState OutS>
class RetaggingBarrierCollectorTask
    : public hh::AbstractTask<1, MeshData<InS>, MeshData<OutS>> {
public:
    template<typename Fn>
    RetaggingBarrierCollectorTask(int nmeshes, std::string name,
                                  std::string routines, Fn&& fn)
        : hh::AbstractTask<1, MeshData<InS>, MeshData<OutS>>(std::move(name), 1),
          nmeshes_(nmeshes), routines_(std::move(routines)),
          fn_(std::forward<Fn>(fn)),
          nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData<InS>> data) override {
        collected_[data->nm - nmOffset_] = data;
        if (++count_ == nmeshes_) {
            auto t0 = std::chrono::steady_clock::now();
            fn_(collected_);
            auto t1 = std::chrono::steady_clock::now();
            totalTime_ += std::chrono::duration<double>(t1 - t0).count();
            ++invocations_;
            count_ = 0;
            for (auto &md : collected_) {
                this->addResult(retag<OutS>(std::move(md)));
            }
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
    BarrierFn<InS> fn_;
    double totalTime_ = 0.0;
    int invocations_ = 0;
    std::vector<std::shared_ptr<MeshData<InS>>> collected_;
};

/// Helper to create a RetaggingBarrierCollectorTask.
template<MeshState InS, MeshState OutS, typename Fn>
inline auto makeRetaggingBarrier(int nmeshes, std::string name,
                                  std::string routines, Fn&& fn) {
    return std::make_shared<RetaggingBarrierCollectorTask<InS, OutS>>(
        nmeshes, std::move(name), std::move(routines), std::forward<Fn>(fn));
}

/// Retagging barrier with TerminationData support for cycle breaking.
/// Same as RetaggingBarrierCollectorTask but also accepts TerminationData
/// and overrides canTerminate() to break structural cycles at shutdown.
template<MeshState InS, MeshState OutS>
class TerminableRetaggingBarrierTask
    : public hh::AbstractTask<2, MeshData<InS>, TerminationData, MeshData<OutS>> {
public:
    template<typename Fn>
    TerminableRetaggingBarrierTask(int nmeshes, std::string name,
                                    std::string routines, Fn&& fn)
        : hh::AbstractTask<2, MeshData<InS>, TerminationData, MeshData<OutS>>(
              std::move(name), 1),
          nmeshes_(nmeshes), routines_(std::move(routines)),
          fn_(std::forward<Fn>(fn)),
          nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData<InS>> data) override {
        collected_[data->nm - nmOffset_] = data;
        if (++count_ == nmeshes_) {
            auto t0 = std::chrono::steady_clock::now();
            fn_(collected_);
            auto t1 = std::chrono::steady_clock::now();
            totalTime_ += std::chrono::duration<double>(t1 - t0).count();
            ++invocations_;
            count_ = 0;
            for (auto &md : collected_) {
                this->addResult(retag<OutS>(std::move(md)));
            }
        }
    }

    void execute(std::shared_ptr<TerminationData>) override { done_ = true; }

    [[nodiscard]] bool canTerminate() const override {
        return done_;
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
    bool done_ = false;
    int nmeshes_, nmOffset_, count_ = 0;
    std::string routines_;
    BarrierFn<InS> fn_;
    double totalTime_ = 0.0;
    int invocations_ = 0;
    std::vector<std::shared_ptr<MeshData<InS>>> collected_;
};

template<MeshState InS, MeshState OutS, typename Fn>
inline auto makeTerminableRetaggingBarrier(int nmeshes, std::string name,
                                            std::string routines, Fn&& fn) {
    return std::make_shared<TerminableRetaggingBarrierTask<InS, OutS>>(
        nmeshes, std::move(name), std::move(routines), std::forward<Fn>(fn));
}

/// Eager dual-input barrier with TerminationData support and retagged output.
/// Same as EagerDualInputBarrierTask but accepts TerminationData for cycle
/// breaking and emits MeshData<OutS> instead of MeshData<>.
template<MeshState OutS>
class TerminableEagerDualInputBarrierTask
    : public hh::AbstractTask<3, MeshData<>, BarrierData, TerminationData, MeshData<OutS>> {
public:
    template<typename Fn>
    TerminableEagerDualInputBarrierTask(int nmeshes, std::string name,
                                         std::string routines, Fn&& fn)
        : hh::AbstractTask<3, MeshData<>, BarrierData, TerminationData, MeshData<OutS>>(
              std::move(name), 1),
          nmeshes_(nmeshes), routines_(std::move(routines)),
          fn_(std::forward<Fn>(fn)),
          nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData<>> data) override {
        collected_[data->nm - nmOffset_] = data;
        if (++meshCount_ == nmeshes_) {
            auto t0 = std::chrono::steady_clock::now();
            fn_(collected_);
            auto t1 = std::chrono::steady_clock::now();
            totalTime_ += std::chrono::duration<double>(t1 - t0).count();
            ++invocations_;
            fnDone_ = true;
            tryEmit();
        }
    }

    void execute(std::shared_ptr<BarrierData>) override {
        exchDone_ = true;
        tryEmit();
    }

    void execute(std::shared_ptr<TerminationData>) override { done_ = true; }

    [[nodiscard]] bool canTerminate() const override {
        return done_;
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
    void tryEmit() {
        if (fnDone_ && exchDone_) {
            fnDone_ = false;
            exchDone_ = false;
            meshCount_ = 0;
            for (auto &md : collected_) {
                this->addResult(retag<OutS>(md));
                md = nullptr;
            }
        }
    }

    bool done_ = false;
    int nmeshes_, nmOffset_, meshCount_ = 0;
    bool fnDone_ = false, exchDone_ = false;
    std::string routines_;
    BarrierFn<> fn_;
    double totalTime_ = 0.0;
    int invocations_ = 0;
    std::vector<std::shared_ptr<MeshData<>>> collected_;
};

template<MeshState OutS, typename Fn>
inline auto makeTerminableEagerDualInputBarrier(int nmeshes, std::string name,
                                                 std::string routines, Fn&& fn) {
    return std::make_shared<TerminableEagerDualInputBarrierTask<OutS>>(
        nmeshes, std::move(name), std::move(routines), std::forward<Fn>(fn));
}

/// Eager dual-mesh barrier with TerminationData support.
/// Collects N MeshData<> (primary) + N MeshData<InS> (secondary).
/// Runs barrier function eagerly when all N primary meshes arrive.
/// Emits N MeshData<OutS> only when BOTH sets are complete.
template<MeshState InS, MeshState OutS>
class TerminableEagerDualMeshBarrierTask
    : public hh::AbstractTask<3, MeshData<>, MeshData<InS>, TerminationData, MeshData<OutS>> {
public:
    template<typename Fn>
    TerminableEagerDualMeshBarrierTask(int nmeshes, std::string name,
                                       std::string routines, Fn&& fn)
        : hh::AbstractTask<3, MeshData<>, MeshData<InS>, TerminationData, MeshData<OutS>>(
              std::move(name), 1),
          nmeshes_(nmeshes), routines_(std::move(routines)),
          fn_(std::forward<Fn>(fn)),
          nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData<>> data) override {
        collected_[data->nm - nmOffset_] = data;
        if (++primaryCount_ == nmeshes_) {
            auto t0 = std::chrono::steady_clock::now();
            fn_(collected_);
            auto t1 = std::chrono::steady_clock::now();
            totalTime_ += std::chrono::duration<double>(t1 - t0).count();
            ++invocations_;
            fnDone_ = true;
            tryEmit();
        }
    }

    void execute(std::shared_ptr<MeshData<InS>>) override {
        if (++secondaryCount_ == nmeshes_) {
            exchDone_ = true;
            tryEmit();
        }
    }

    void execute(std::shared_ptr<TerminationData>) override { done_ = true; }

    [[nodiscard]] bool canTerminate() const override { return done_; }

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
    void tryEmit() {
        if (fnDone_ && exchDone_) {
            fnDone_ = false;
            exchDone_ = false;
            primaryCount_ = 0;
            secondaryCount_ = 0;
            for (auto &md : collected_) {
                this->addResult(retag<OutS>(md));
                md = nullptr;
            }
        }
    }

    bool done_ = false;
    int nmeshes_, nmOffset_;
    int primaryCount_ = 0, secondaryCount_ = 0;
    bool fnDone_ = false, exchDone_ = false;
    std::string routines_;
    BarrierFn<> fn_;
    double totalTime_ = 0.0;
    int invocations_ = 0;
    std::vector<std::shared_ptr<MeshData<>>> collected_;
};

template<MeshState InS, MeshState OutS, typename Fn>
inline auto makeTerminableEagerDualMeshBarrier(int nmeshes, std::string name,
                                                std::string routines, Fn&& fn) {
    return std::make_shared<TerminableEagerDualMeshBarrierTask<InS, OutS>>(
        nmeshes, std::move(name), std::move(routines), std::forward<Fn>(fn));
}

#endif // BARRIER_STATE_H
