#ifndef CHANGE_TIMESTEP_STATE_H
#define CHANGE_TIMESTEP_STATE_H

#include <hedgehog/hedgehog.h>
#include "../data/change_timestep_data.h"
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/termination_data.h"
#include "../fds_fortran_interface.h"

/// Merged collector + retry loop state.
///
/// Collects N MeshData<> tokens from the parallel kernel (or receives
/// RetrySequenceData from the bypass path), runs post-kernel work,
/// checks for CFL retry, and either cycles back or exits.
///
/// Output types (Hedgehog type-based routing):
///   - RetrySequenceData → cycles back to RetryPreKernel for another retry
///   - BarrierData → exits the subgraph (retry complete or no retry needed)
class RetryLoopState : public hh::AbstractState<3, MeshData<>, RetrySequenceData, TerminationData,
                                                  RetrySequenceData, BarrierData> {
public:
    explicit RetryLoopState(int nmeshes)
        : nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    /// Collect MeshData<> from parallel kernel (N tokens).
    void execute(std::shared_ptr<MeshData<>> data) override {
        collected_[data->nm - nmOffset_] = data;
        if (++count_ == nmeshes_) {
            count_ = 0;
            processPostKernel();
        }
    }

    /// Bypass from RetryPreKernel (done=true → no retry needed).
    void execute(std::shared_ptr<RetrySequenceData> data) override {
        if (data->done) {
            auto bd = std::make_shared<BarrierData>();
            bd->meshes = data->meshes;
            this->addResult(bd);
        }
    }

    /// Handle termination signal from graph input.
    void execute(std::shared_ptr<TerminationData>) override {
        done_ = true;
    }

    [[nodiscard]] bool isDone() const { return done_; }

private:
    void processPostKernel() {
        double t = collected_[0]->t;
        double dt = collected_[0]->dt;

        // Post-kernel work (was RetryPostKernelTask)
        fds_exchange_divergence_info();

        for (auto &md : collected_) {
            fds_divergence_part_2(dt, md->nm);
        }

        fds_pressure_iteration(t, dt);
        fds_init_change_time_step(dt);

        for (auto &md : collected_) {
            fds_velocity_predictor(t + dt, dt, md->nm);
        }

        fds_stop_check_zero();

        // Check if we should stop due to instability
        int stopStatus = fds_get_stop_status();
        if (stopStatus != 0) {
            emitExit();
            return;
        }

        // Check if another retry is needed
        int needRetry = 0;
        double newDt = 0.0;
        fds_check_change_time_step(&needRetry, &newDt);

        if (!needRetry) {
            emitExit();
        } else {
            ++iteration_;
            auto retryData = std::make_shared<RetrySequenceData>(
                collected_, t, newDt, iteration_, false);
            this->addResult(retryData);
        }
    }

    void emitExit() {
        auto bd = std::make_shared<BarrierData>();
        bd->meshes.reserve(nmeshes_);
        for (auto &md : collected_) {
            bd->meshes.push_back(md);
            md = nullptr;
        }
        iteration_ = 0;
        this->addResult(bd);
    }

    int nmeshes_, nmOffset_, count_ = 0, iteration_ = 0;
    bool done_ = false;
    std::vector<std::shared_ptr<MeshData<>>> collected_;
};

/// Custom state manager for the retry loop cycle.
class RetryLoopStateManager
    : public hh::StateManager<3, MeshData<>, RetrySequenceData, TerminationData,
                               RetrySequenceData, BarrierData> {
public:
    RetryLoopStateManager(
        std::shared_ptr<RetryLoopState> const &state,
        std::string const &name)
        : hh::StateManager<3, MeshData<>, RetrySequenceData, TerminationData,
                            RetrySequenceData, BarrierData>(
              state, name) {}

    [[nodiscard]] bool canTerminate() const override {
        this->state()->lock();
        auto s = std::dynamic_pointer_cast<RetryLoopState>(this->state());
        bool ret = s->isDone();
        this->state()->unlock();
        return ret;
    }
};

#endif // CHANGE_TIMESTEP_STATE_H
