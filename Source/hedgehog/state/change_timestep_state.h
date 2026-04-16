#ifndef CHANGE_TIMESTEP_STATE_H
#define CHANGE_TIMESTEP_STATE_H

#include <hedgehog/hedgehog.h>
#include "../data/change_timestep_data.h"
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/termination_data.h"
#include "../fds_fortran_interface.h"

/// Retry check state: collects VelPred kernel results, checks CFL retry.
///
/// Three input paths:
///   1. MeshData<> from VelocityPredictor kernel (N tokens per iteration).
///      Collects N, runs STOP_CHECK + retry check.
///   2. RetrySequenceData(done=true) from RetryPreKernel bypass (no retry needed).
///      Runs CC_END_STEP + MESH_EXCHANGE(3), emits MeshData<>.
///   3. TerminationData for clean shutdown.
///
/// Output types (Hedgehog type-based routing):
///   - RetrySequenceData → cycles back to RetryPreKernel for another retry
///   - MeshData<> → exits the subgraph (retry complete or no retry needed)
class RetryCheckState : public hh::AbstractState<3, MeshData<>, RetrySequenceData, TerminationData,
                                                  RetrySequenceData, MeshData<>> {
public:
    RetryCheckState(int nmeshes, bool ccIBM)
        : nmeshes_(nmeshes), ccIBM_(ccIBM),
          nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    /// Collect MeshData<> from VelocityPredictor kernel (N tokens).
    void execute(std::shared_ptr<MeshData<>> data) override {
        collected_[data->nm - nmOffset_] = data;
        if (++count_ == nmeshes_) {
            count_ = 0;
            processCheck();
        }
    }

    /// Bypass from RetryPreKernel (done=true → no retry needed).
    void execute(std::shared_ptr<RetrySequenceData> data) override {
        if (data->done) {
            if (ccIBM_) { fds_cc_end_step(data->meshes[0]->t, data->meshes[0]->dt, 0); }
            fds_mesh_exchange(3);
            for (auto &md : data->meshes) {
                this->addResult(md);
            }
        }
    }

    /// Handle termination signal from graph input.
    void execute(std::shared_ptr<TerminationData>) override {
        done_ = true;
    }

    [[nodiscard]] bool isDone() const { return done_; }

private:
    void processCheck() {
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
            std::vector<std::shared_ptr<MeshData<>>> meshes(collected_.begin(), collected_.end());
            std::fill(collected_.begin(), collected_.end(), nullptr);
            auto retryData = std::make_shared<RetrySequenceData>(
                meshes, meshes[0]->t, newDt, iteration_, false);
            this->addResult(retryData);
        }
    }

    void emitExit() {
        if (ccIBM_) { fds_cc_end_step(collected_[0]->t, collected_[0]->dt, 0); }
        fds_mesh_exchange(3);
        for (auto &md : collected_) {
            this->addResult(md);
            md = nullptr;
        }
        iteration_ = 0;
    }

    int nmeshes_, nmOffset_, count_ = 0, iteration_ = 0;
    bool ccIBM_;
    bool done_ = false;
    std::vector<std::shared_ptr<MeshData<>>> collected_;
};

/// Custom state manager for the retry check cycle.
class RetryCheckStateManager
    : public hh::StateManager<3, MeshData<>, RetrySequenceData, TerminationData,
                               RetrySequenceData, MeshData<>> {
public:
    RetryCheckStateManager(
        std::shared_ptr<RetryCheckState> const &state,
        std::string const &name)
        : hh::StateManager<3, MeshData<>, RetrySequenceData, TerminationData,
                            RetrySequenceData, MeshData<>>(
              state, name) {}

    [[nodiscard]] bool canTerminate() const override {
        this->state()->lock();
        auto s = std::dynamic_pointer_cast<RetryCheckState>(this->state());
        bool ret = s->isDone();
        this->state()->unlock();
        return ret;
    }
};

#endif // CHANGE_TIMESTEP_STATE_H
