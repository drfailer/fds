#ifndef CHANGE_TIMESTEP_STATE_H
#define CHANGE_TIMESTEP_STATE_H

#include <hedgehog/hedgehog.h>
#include "../data/change_timestep_data.h"
#include "../data/mesh_data.h"
#include "../data/termination_data.h"
#include "../fds_fortran_interface.h"

/// State that manages the retry loop cycle.
///
/// Receives RetrySequenceData after each retry sequence pass.
/// Includes the post-kernel work (divergence exchange, div_p2, pressure,
/// velocity predictor) that was previously in RetryPostKernelTask.
///
/// Two output types for Hedgehog type-based routing:
///   - RetrySequenceData → cycles back to RetryPreKernel for another retry
///   - MeshData → exits the subgraph (retry complete or no retry needed)
class RetryLoopState : public hh::AbstractState<2, RetrySequenceData, TerminationData, RetrySequenceData, MeshData> {
public:
    RetryLoopState() = default;

    void execute(std::shared_ptr<RetrySequenceData> data) override {
        if (data->done) {
            // No retry needed: emit MeshData to exit
            for (auto &md : data->meshes) {
                this->addResult(md);
            }
            return;
        }

        // Post-kernel work (was RetryPostKernelTask)
        fds_exchange_divergence_info();

        for (auto &md : data->meshes) {
            fds_divergence_part_2(data->dt, md->nm);
        }

        fds_pressure_iteration(data->t, data->dt);
        fds_init_change_time_step(data->dt);

        for (auto &md : data->meshes) {
            fds_velocity_predictor(data->t + data->dt, data->dt, md->nm);
        }

        fds_stop_check_zero();

        // Check if we should stop due to instability
        int stopStatus = fds_get_stop_status();
        if (stopStatus != 0) {
            for (auto &md : data->meshes) {
                this->addResult(md);
            }
            return;
        }

        // Check if another retry is needed
        int needRetry = 0;
        double newDt = 0.0;
        fds_check_change_time_step(&needRetry, &newDt);

        if (!needRetry) {
            for (auto &md : data->meshes) {
                this->addResult(md);
            }
        } else {
            auto retryData = std::make_shared<RetrySequenceData>(
                data->meshes, data->t, newDt, data->iteration + 1, false);
            this->addResult(retryData);
        }
    }

    /// Handle termination signal from graph input.
    void execute(std::shared_ptr<TerminationData>) override {
        done_ = true;
    }

    [[nodiscard]] bool isDone() const { return done_; }

private:
    bool done_ = false;
};

/// Custom state manager for the retry loop cycle.
class RetryLoopStateManager
    : public hh::StateManager<2, RetrySequenceData, TerminationData, RetrySequenceData, MeshData> {
public:
    RetryLoopStateManager(
        std::shared_ptr<RetryLoopState> const &state,
        std::string const &name)
        : hh::StateManager<2, RetrySequenceData, TerminationData, RetrySequenceData, MeshData>(
              state, name) {}

    [[nodiscard]] bool canTerminate() const override {
        this->state()->lock();
        auto s = std::dynamic_pointer_cast<RetryLoopState>(this->state());
        bool ret = s->isDone();
        this->state()->unlock();
        return ret;
    }
};

/// Collector that gathers N MeshData tokens after parallel retry kernel
/// processing and reconstructs a RetrySequenceData for the pipeline.
class RetryMomentumDivCollector
    : public hh::AbstractState<1, MeshData, RetrySequenceData> {
public:
    explicit RetryMomentumDivCollector(int nmeshes)
        : nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_[data->nm - nmOffset_] = data;
        ++count_;

        if (count_ == nmeshes_) {
            auto retryData = std::make_shared<RetrySequenceData>(
                collected_, collected_[0]->t, collected_[0]->dt, 0, false);
            collected_.resize(nmeshes_, nullptr);
            count_ = 0;
            this->addResult(retryData);
        }
    }

private:
    int nmeshes_;
    int nmOffset_;
    int count_ = 0;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

#endif // CHANGE_TIMESTEP_STATE_H
