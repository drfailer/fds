#ifndef PRESSURE_CONVERGENCE_STATE_H
#define PRESSURE_CONVERGENCE_STATE_H

#include <hedgehog/hedgehog.h>
#include <chrono>
#include <iomanip>
#include <sstream>
#include "../data/pressure_iteration_data.h"
#include "../data/mesh_data.h"
#include "../data/termination_signal.h"
#include "../fds_fortran_interface.h"

/// Pressure iteration convergence barrier.
///
/// Collects N MeshData tokens from PressureSolveKernel, then runs:
///   1. fds_mesh_exchange(5) — exchange FVX/FVY/FVZ/H between neighbors
///   2. velocity_error_kernel — compute velocity error per mesh
///   3. convergence_check — MPI reduction + tolerance check
///
/// Routes:
///   - Converged: emit N MeshData tokens (exit subgraph)
///   - Not converged: call increment(), emit PressureIterMeshData (cycle back)
///
/// Also handles ITERATE_PRESSURE=false (always exit after first pass).
class PressureConvergenceState
    : public hh::AbstractState<1, MeshData, PressureIterMeshData, MeshData> {
public:
    PressureConvergenceState(int nmeshes, double tEnd, bool predictor,
                             std::shared_ptr<TerminationSignal> termSignal)
        : nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()),
          tEnd_(tEnd), predictor_(predictor),
          termSignal_(std::move(termSignal)) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_[data->nm - nmOffset_] = data;

        if (++count_ == nmeshes_) {
            double t = collected_[0]->t;
            double dt = collected_[0]->dt;
            lastT_ = t;
            lastDt_ = dt;

            int converged;
            if (fds_iterate_pressure()) {
                auto t0 = std::chrono::steady_clock::now();
                fds_mesh_exchange(5);
                for (int i = 0; i < nmeshes_; ++i) {
                    fds_compute_velocity_error_kernel(collected_[i]->nm, dt);
                }
                fds_pressure_iteration_check_convergence(t, dt);
                auto t1 = std::chrono::steady_clock::now();
                convTime_ += std::chrono::duration<double>(t1 - t0).count();
                converged = fds_pressure_iteration_converged();
            } else {
                // ITERATE_PRESSURE is false: always exit after one pass
                converged = 1;
            }

            ++invocations_;
            lastConverged_ = (converged != 0);

            if (converged) {
                if (predictor_) {
                    fds_init_change_time_step(dt);
                }
                for (int i = 0; i < nmeshes_; ++i) {
                    this->addResult(std::move(collected_[i]));
                }
            } else {
                // Increment counter for next iteration
                fds_pressure_iteration_increment();

                for (int i = 0; i < nmeshes_; ++i) {
                    this->addResult(std::make_shared<PressureIterMeshData>(std::move(collected_[i])));
                }
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
        oss << "MESH_EXCHANGE(5)\\nVEL_ERROR\\nCONVERGENCE_CHECK\\n"
            << std::fixed << std::setprecision(3)
            << "conv " << convTime_ << "s"
            << " / " << invocations_ << " calls";
        if (invocations_ > 0) {
            oss << " / avg " << std::setprecision(3)
                << (convTime_ * 1000.0 / invocations_) << "ms";
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
    double convTime_ = 0.0;
    int invocations_ = 0;
};

/// Custom state manager for PressureConvergenceState with canTerminate.
class PressureConvergenceManager
    : public hh::StateManager<1, MeshData, PressureIterMeshData, MeshData> {
public:
    PressureConvergenceManager(
        std::shared_ptr<PressureConvergenceState> const& state,
        std::string const& name)
        : hh::StateManager<1, MeshData, PressureIterMeshData, MeshData>(
              state, name) {}

    [[nodiscard]] bool canTerminate() const override {
        this->state()->lock();
        auto s = std::dynamic_pointer_cast<PressureConvergenceState>(
            this->state());
        bool ret = (s->reachedEnd() && s->lastConverged()) || s->isTerminated();
        this->state()->unlock();
        return ret;
    }

    [[nodiscard]] std::string extraPrintingInformation() const override {
        this->state()->lock();
        auto ret = std::dynamic_pointer_cast<PressureConvergenceState>(
            this->state())->info();
        this->state()->unlock();
        return ret;
    }
};

#endif // PRESSURE_CONVERGENCE_STATE_H
