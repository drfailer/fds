#ifndef PRESSURE_CONVERGENCE_STATE_H
#define PRESSURE_CONVERGENCE_STATE_H

#include <hedgehog/hedgehog.h>
#include <chrono>
#include <iomanip>
#include <sstream>
#include "../data/pressure_iteration_data.h"
#include "../data/mesh_data.h"
#include "../data/termination_data.h"
#include "../fds_fortran_interface.h"

/// Pressure iteration convergence check barrier.
///
/// Collects N MeshData tokens (after post-solve exchange and velocity error),
/// then runs the MPI convergence check.
///
/// Routes:
///   - Converged: emit N MeshData tokens (exit subgraph)
///   - Not converged: call increment(), emit PressureIterMeshData (cycle back)
///
/// Also handles ITERATE_PRESSURE=false (always exit after first pass).
///
/// Termination: receives TerminationData from the graph input when the
/// simulation is complete, setting done_=true so canTerminate() returns true.
class PressureConvergenceState
    : public hh::AbstractState<2, MeshData, TerminationData, PressureIterMeshData, MeshData> {
public:
    /// @param nmeshes  Number of local meshes
    /// @param predictor True for predictor phase
    PressureConvergenceState(int nmeshes, bool predictor)
        : nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()),
          predictor_(predictor) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_[data->nm - nmOffset_] = data;

        if (++count_ == nmeshes_) {
            double t = collected_[0]->t;
            double dt = collected_[0]->dt;

            int converged;
            if (fds_iterate_pressure()) {
                auto t0 = std::chrono::steady_clock::now();
                fds_pressure_iteration_check_convergence(t, dt);
                auto t1 = std::chrono::steady_clock::now();
                convTime_ += std::chrono::duration<double>(t1 - t0).count();
                converged = fds_pressure_iteration_converged();
            } else {
                // ITERATE_PRESSURE is false: always exit after one pass
                converged = 1;
            }

            ++invocations_;

            if (converged) {
                if (predictor_) {
                    fds_init_change_time_step(dt);
                }
                this->batchAddResult(collected_);
            } else {
                // Increment counter for next iteration
                fds_pressure_iteration_increment();

                for (int i = 0; i < nmeshes_; ++i) {
                    this->bufferResult(std::make_shared<PressureIterMeshData>(std::move(collected_[i])));
                }
                this->flushResults<PressureIterMeshData>();
            }

            count_ = 0;
        }
    }

    /// Handle termination signal from graph input.
    void execute(std::shared_ptr<TerminationData>) override {
        done_ = true;
    }

    [[nodiscard]] bool isDone() const { return done_; }

    [[nodiscard]] std::string info() const {
        std::ostringstream oss;
        oss << "CONVERGENCE_CHECK\\n"
            << std::fixed << std::setprecision(3)
            << "check " << convTime_ << "s"
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
    bool predictor_;
    bool done_ = false;
    double convTime_ = 0.0;
    int invocations_ = 0;
};

/// Custom state manager for PressureConvergenceState with canTerminate.
class PressureConvergenceManager
    : public hh::StateManager<2, MeshData, TerminationData, PressureIterMeshData, MeshData> {
public:
    PressureConvergenceManager(
        std::shared_ptr<PressureConvergenceState> const& state,
        std::string const& name)
        : hh::StateManager<2, MeshData, TerminationData, PressureIterMeshData, MeshData>(
              state, name) {}

    [[nodiscard]] bool canTerminate() const override {
        this->state()->lock();
        auto s = std::dynamic_pointer_cast<PressureConvergenceState>(
            this->state());
        bool ret = s->isDone();
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
