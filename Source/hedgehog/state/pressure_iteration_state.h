#ifndef PRESSURE_ITERATION_STATE_H
#define PRESSURE_ITERATION_STATE_H

#include <hedgehog/hedgehog.h>
#include "../data/pressure_iteration_data.h"
#include "../data/mesh_data.h"
#include "../data/termination_signal.h"
#include "../fds_fortran_interface.h"

/// Merged collector + post-loop state for pressure iteration.
///
/// Collects N MeshData tokens from the parallel pressure solve kernel,
/// runs Phase 3 sequential work (MESH_EXCHANGE(5) + velocity error +
/// convergence check), then routes:
///   - PressureIterData -> cycles back to PressurePreKernelTask (not converged)
///   - MeshData -> exits the sub-graph (converged)
///
/// This replaces the former PressureSolveCollector + PressurePostLoopState
/// two-node chain, eliminating one queue transition per pressure iteration.
class PressurePostCollector
    : public hh::AbstractState<1, MeshData, PressureIterData, MeshData> {
public:
    PressurePostCollector(int nmeshes, double tEnd, bool predictor,
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
                fds_mesh_exchange(5);
                for (auto& md : collected_) {
                    fds_compute_velocity_error_kernel(md->nm, dt);
                }
                fds_pressure_iteration_check_convergence(t, dt);
            }

            int converged = fds_pressure_iteration_converged();
            lastConverged_ = (converged != 0);

            if (converged) {
                if (predictor_) {
                    fds_init_change_time_step(dt);
                }
                for (auto& md : collected_) {
                    this->addResult(md);
                    md = nullptr;
                }
            } else {
                // Build PressureIterData for cycle
                std::vector<std::shared_ptr<MeshData>> meshes;
                meshes.reserve(nmeshes_);
                for (auto& md : collected_) {
                    meshes.push_back(md);
                    md = nullptr;
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
};

/// Custom state manager for the merged pressure post-collector.
class PressurePostCollectorManager
    : public hh::StateManager<1, MeshData, PressureIterData, MeshData> {
public:
    PressurePostCollectorManager(
        std::shared_ptr<PressurePostCollector> const& state,
        std::string const& name)
        : hh::StateManager<1, MeshData, PressureIterData, MeshData>(
              state, name) {}

    [[nodiscard]] bool canTerminate() const override {
        this->state()->lock();
        auto s = std::dynamic_pointer_cast<PressurePostCollector>(
            this->state());
        bool ret = (s->reachedEnd() && s->lastConverged()) || s->isTerminated();
        this->state()->unlock();
        return ret;
    }
};

#endif // PRESSURE_ITERATION_STATE_H
