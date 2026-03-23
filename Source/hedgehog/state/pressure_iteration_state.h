#ifndef PRESSURE_ITERATION_STATE_H
#define PRESSURE_ITERATION_STATE_H

#include <hedgehog/hedgehog.h>
#include <chrono>
#include <iomanip>
#include <sstream>
#include "../data/pressure_iteration_data.h"
#include "../data/mesh_data.h"
#include "../data/termination_signal.h"
#include "../fds_fortran_interface.h"

/// Merged collector + pre-kernel state for pressure iteration.
///
/// Collects N MeshData tokens from upstream (initial entry) or unpacks
/// PressureIterData (cycle), runs Phase 1 sequential work (pressure iteration
/// init + baroclinic correction + mesh exchange), then scatters MeshData to
/// the parallel solve kernel.
///
/// Replaces the former PredPressureCollector (CollectorState) +
/// PressurePreKernelTask two-node chain, eliminating one queue transition
/// and the BarrierData intermediate.
class PressurePreCollector
    : public hh::AbstractState<2, MeshData, PressureIterData, MeshData> {
public:
    explicit PressurePreCollector(int nmeshes)
        : nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    /// Initial entry: collect N MeshData tokens, then run init + Phase 1.
    void execute(std::shared_ptr<MeshData> data) override {
        collected_[data->nm - nmOffset_] = data;
        if (++count_ == nmeshes_) {
            auto t0 = std::chrono::steady_clock::now();
            fds_pressure_iteration_init();
            auto t1 = std::chrono::steady_clock::now();
            initTime_ += std::chrono::duration<double>(t1 - t0).count();

            runPhase1();
        }
    }

    /// Cycle entry: unpack PressureIterData, run Phase 1 (no init).
    void execute(std::shared_ptr<PressureIterData> data) override {
        for (auto& md : data->meshes) {
            collected_[md->nm - nmOffset_] = md;
        }
        runPhase1();
    }

    [[nodiscard]] std::string info() const {
        std::ostringstream oss;
        oss << "PRESSURE_ITERATION_INIT\\n"
            << "BAROCLINIC_CORRECTION\\n"
            << "MESH_EXCHANGE(5)\\n"
            << std::fixed << std::setprecision(3)
            << "init " << initTime_ << "s"
            << " / baro " << baroTime_ << "s"
            << " / exch " << exchTime_ << "s\\n"
            << invocations_ << " calls";
        if (invocations_ > 0) {
            double total = initTime_ + baroTime_ + exchTime_;
            oss << " / avg " << std::setprecision(3)
                << (total * 1000.0 / invocations_) << "ms";
        }
        return oss.str();
    }

private:
    void runPhase1() {
        double t = collected_[0]->t;

        fds_pressure_iteration_increment();

        if (fds_pressure_iteration_needs_baroclinic()) {
            auto t0 = std::chrono::steady_clock::now();
            for (int i = 0; i < nmeshes_; ++i) {
                fds_baroclinic_correction(t, collected_[i]->nm);
            }
            auto t1 = std::chrono::steady_clock::now();
            baroTime_ += std::chrono::duration<double>(t1 - t0).count();

            fds_mesh_exchange(5);
            auto t2 = std::chrono::steady_clock::now();
            exchTime_ += std::chrono::duration<double>(t2 - t1).count();
        }

        ++invocations_;
        count_ = 0;

        for (int i = 0; i < nmeshes_; ++i) {
            this->addResult(collected_[i]);
            collected_[i] = nullptr;
        }
    }

    int nmeshes_, nmOffset_, count_ = 0;
    std::vector<std::shared_ptr<MeshData>> collected_;
    double initTime_ = 0.0, baroTime_ = 0.0, exchTime_ = 0.0;
    int invocations_ = 0;
};

/// StateManager for PressurePreCollector with profiling.
class PressurePreCollectorManager
    : public hh::StateManager<2, MeshData, PressureIterData, MeshData> {
public:
    PressurePreCollectorManager(
        std::shared_ptr<PressurePreCollector> const& state,
        std::string const& name)
        : hh::StateManager<2, MeshData, PressureIterData, MeshData>(
              state, name) {}

    [[nodiscard]] std::string extraPrintingInformation() const override {
        this->state()->lock();
        auto ret = std::dynamic_pointer_cast<PressurePreCollector>(
            this->state())->info();
        this->state()->unlock();
        return ret;
    }
};


/// Merged collector + post-loop state for pressure iteration.
///
/// Collects N MeshData tokens from the parallel pressure solve kernel,
/// runs Phase 3 sequential work (MESH_EXCHANGE(5) + velocity error +
/// convergence check), then routes:
///   - PressureIterData -> cycles back to PressurePreCollector (not converged)
///   - MeshData -> exits the sub-graph (converged)
class PressurePostLoop
    : public hh::AbstractState<1, MeshData, PressureIterData, MeshData> {
public:
    PressurePostLoop(int nmeshes, double tEnd, bool predictor,
                     std::shared_ptr<TerminationSignal> termSignal)
        : nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()),
          tEnd_(tEnd), predictor_(predictor),
          termSignal_(std::move(termSignal)) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_[data->nm - nmOffset_] = data;
        ++count_;

        if (count_ == nmeshes_) {
            double t = collected_[0]->t;
            double dt = collected_[0]->dt;
            lastT_ = t;
            lastDt_ = dt;

            // Phase 3: sequential post-kernel work
            int iteratePressure = fds_iterate_pressure();
            if (iteratePressure) {
                auto t0 = std::chrono::steady_clock::now();
                fds_mesh_exchange(5);
                auto t1 = std::chrono::steady_clock::now();
                exchTime_ += std::chrono::duration<double>(t1 - t0).count();

                for (int i = 0; i < nmeshes_; ++i) {
                    fds_compute_velocity_error_kernel(collected_[i]->nm, dt);
                }
                auto t2 = std::chrono::steady_clock::now();
                velErrTime_ += std::chrono::duration<double>(t2 - t1).count();

                fds_pressure_iteration_check_convergence(t, dt);
                auto t3 = std::chrono::steady_clock::now();
                convTime_ += std::chrono::duration<double>(t3 - t2).count();
            }

            ++invocations_;
            int converged = fds_pressure_iteration_converged();
            lastConverged_ = (converged != 0);

            if (converged) {
                if (predictor_) {
                    fds_init_change_time_step(dt);
                }
                for (int i = 0; i < nmeshes_; ++i) {
                    this->addResult(collected_[i]);
                    collected_[i] = nullptr;
                }
            } else {
                std::vector<std::shared_ptr<MeshData>> meshes;
                meshes.reserve(nmeshes_);
                for (int i = 0; i < nmeshes_; ++i) {
                    meshes.push_back(collected_[i]);
                    collected_[i] = nullptr;
                }
                this->addResult(std::make_shared<PressureIterData>(
                    std::move(meshes), t, dt));
            }

            count_ = 0;
        }
    }

    [[nodiscard]] bool reachedEnd() const {
        if (predictor_) {
            return lastT_ + lastDt_ >= tEnd_;
        } else {
            return lastT_ >= tEnd_;
        }
    }

    [[nodiscard]] bool isTerminated() const {
        return termSignal_->isTerminated();
    }

    [[nodiscard]] bool lastConverged() const { return lastConverged_; }

    [[nodiscard]] std::string info() const {
        std::ostringstream oss;
        oss << "MESH_EXCHANGE(5)\\n"
            << "COMPUTE_VELOCITY_ERROR\\n"
            << "CONVERGENCE_CHECK\\n"
            << std::fixed << std::setprecision(3)
            << "exch " << exchTime_ << "s"
            << " / vel_err " << velErrTime_ << "s"
            << " / conv " << convTime_ << "s\\n"
            << invocations_ << " calls";
        if (invocations_ > 0) {
            double total = exchTime_ + velErrTime_ + convTime_;
            oss << " / avg " << std::setprecision(3)
                << (total * 1000.0 / invocations_) << "ms";
        }
        return oss.str();
    }

private:
    int nmeshes_;
    int nmOffset_;
    int count_ = 0;
    std::vector<std::shared_ptr<MeshData>> collected_;
    double tEnd_;
    bool predictor_;
    std::shared_ptr<TerminationSignal> termSignal_;
    double lastT_ = 0.0;
    double lastDt_ = 0.0;
    bool lastConverged_ = false;
    double exchTime_ = 0.0, velErrTime_ = 0.0, convTime_ = 0.0;
    int invocations_ = 0;
};

/// Custom state manager for PressurePostLoop with canTerminate and profiling.
class PressurePostLoopManager
    : public hh::StateManager<1, MeshData, PressureIterData, MeshData> {
public:
    PressurePostLoopManager(
        std::shared_ptr<PressurePostLoop> const& state,
        std::string const& name)
        : hh::StateManager<1, MeshData, PressureIterData, MeshData>(
              state, name) {}

    [[nodiscard]] bool canTerminate() const override {
        this->state()->lock();
        auto s = std::dynamic_pointer_cast<PressurePostLoop>(
            this->state());
        bool ret = (s->reachedEnd() && s->lastConverged()) || s->isTerminated();
        this->state()->unlock();
        return ret;
    }

    [[nodiscard]] std::string extraPrintingInformation() const override {
        this->state()->lock();
        auto ret = std::dynamic_pointer_cast<PressurePostLoop>(
            this->state())->info();
        this->state()->unlock();
        return ret;
    }
};

#endif // PRESSURE_ITERATION_STATE_H
