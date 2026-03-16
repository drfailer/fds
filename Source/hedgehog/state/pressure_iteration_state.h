#ifndef PRESSURE_ITERATION_STATE_H
#define PRESSURE_ITERATION_STATE_H

#include <hedgehog/hedgehog.h>
#include "../data/pressure_iteration_data.h"
#include "../data/mesh_data.h"
#include "../data/termination_signal.h"
#include "../fds_fortran_interface.h"

/// Collector that gathers N MeshData tokens after parallel pressure solve
/// and reconstructs a PressureIterData for the post-kernel pipeline.
class PressureSolveCollector
    : public hh::AbstractState<1, MeshData, PressureIterData> {
public:
    explicit PressureSolveCollector(int nmeshes)
        : nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_[data->nm - nmOffset_] = data;
        ++count_;

        if (count_ == nmeshes_) {
            auto iterData = std::make_shared<PressureIterData>(
                collected_, collected_[0]->t, collected_[0]->dt);
            collected_.assign(nmeshes_, nullptr);
            count_ = 0;
            this->addResult(iterData);
        }
    }

private:
    int nmeshes_;
    int nmOffset_;
    int count_ = 0;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

/// Combined post-kernel + loop state for pressure iteration.
///
/// Receives PressureIterData from the collector, runs Phase 3 sequential work
/// (MESH_EXCHANGE(5) + velocity error + convergence check), then routes:
///   - PressureIterData -> cycles back to PressurePreKernelTask (not converged)
///   - MeshData -> exits the sub-graph (converged)
///
/// Termination uses two conditions checked in canTerminate():
///   1. reachedEnd() && lastConverged(): Data-driven, fires after the last
///      pressure iteration converges on the last timestep. This is the primary
///      mechanism — it triggers during execute(), when Hedgehog re-checks
///      canTerminate(). Without lastConverged(), reachedEnd() alone would fire
///      mid-iteration and terminate the cycle prematurely.
///   2. isTerminated(): Fallback via shared TerminationSignal from the main
///      timestep loop.
class PressurePostLoopState
    : public hh::AbstractState<1, PressureIterData, PressureIterData, MeshData> {
public:
    PressurePostLoopState(double tEnd, bool predictor,
                          std::shared_ptr<TerminationSignal> termSignal)
        : tEnd_(tEnd), predictor_(predictor),
          termSignal_(std::move(termSignal)) {}

    void execute(std::shared_ptr<PressureIterData> data) override {
        lastT_ = data->t;
        lastDt_ = data->dt;

        // Phase 3: sequential post-kernel work
        int iteratePressure = fds_iterate_pressure();
        if (iteratePressure) {
            fds_mesh_exchange(5);
            for (auto& md : data->meshes) {
                fds_compute_velocity_error_kernel(md->nm, data->dt);
            }
            fds_pressure_iteration_check_convergence(data->t, data->dt);
        }

        int converged = fds_pressure_iteration_converged();
        lastConverged_ = (converged != 0);

        if (converged) {
            if (predictor_) {
                fds_init_change_time_step(data->dt);
            }
            for (auto& md : data->meshes) {
                this->addResult(md);
            }
        } else {
            this->addResult(data);
        }
    }

    /// Data-driven termination check.
    /// Predictor: lastT + lastDt >= tEnd (t not yet advanced by PhaseTransition)
    /// Corrector: lastT >= tEnd (t already advanced by PhaseTransition)
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

private:
    double tEnd_;
    bool predictor_;
    std::shared_ptr<TerminationSignal> termSignal_;
    double lastT_ = 0.0;
    double lastDt_ = 0.0;
    bool lastConverged_ = false;
};

/// Custom state manager for the pressure post-loop state.
///
/// canTerminate() requires BOTH reachedEnd() AND lastConverged() to prevent
/// premature termination mid-pressure-iteration on the last timestep.
/// The TerminationSignal serves as a fallback.
class PressurePostLoopStateManager
    : public hh::StateManager<1, PressureIterData, PressureIterData, MeshData> {
public:
    PressurePostLoopStateManager(
        std::shared_ptr<PressurePostLoopState> const& state,
        std::string const& name)
        : hh::StateManager<1, PressureIterData, PressureIterData, MeshData>(
              state, name) {}

    [[nodiscard]] bool canTerminate() const override {
        this->state()->lock();
        auto s = std::dynamic_pointer_cast<PressurePostLoopState>(
            this->state());
        // Primary: data-driven time check AND convergence (prevents mid-iteration termination)
        // Fallback: external termination signal from main timestep loop
        bool ret = (s->reachedEnd() && s->lastConverged()) || s->isTerminated();
        this->state()->unlock();
        return ret;
    }
};

#endif // PRESSURE_ITERATION_STATE_H
