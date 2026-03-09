#ifndef CHANGE_TIMESTEP_STATE_H
#define CHANGE_TIMESTEP_STATE_H

#include <hedgehog/hedgehog.h>
#include "../data/change_timestep_data.h"
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// State that manages the retry loop.
/// After each retry sequence completion, checks if another retry is needed.
/// Either cycles back (done=false) or marks as complete (done=true).
///
/// done_ starts as true because the cycle has not been entered yet.
/// When the cycle is entered (retry needed), done_ is set to false.
/// When the cycle completes (no more retries), done_ is set back to true.
/// This allows canTerminate() to return true in both the "no retry" case
/// (cycle never entered) and the "retry completed" case.
class RetryLoopState : public hh::AbstractState<1, RetrySequenceData, RetrySequenceData> {
public:
    RetryLoopState() = default;

    void execute(std::shared_ptr<RetrySequenceData> data) override {
        // If already marked done (bypass from CheckRetryTask), pass through
        if (data->done) {
            done_ = true;
            this->addResult(data);
            return;
        }

        // Cycle is active
        done_ = false;

        // Check stop status
        int stopStatus = fds_get_stop_status();
        if (stopStatus != 0) {
            // Stop condition met: mark as done and emit
            data->done = true;
            done_ = true;
            this->addResult(data);
            return;
        }

        // Check if another retry is needed
        int needRetry = 0;
        double newDt = 0.0;
        fds_check_change_time_step(&needRetry, &newDt);

        if (!needRetry) {
            // No more retries needed: mark as done and emit
            data->done = true;
            done_ = true;
            this->addResult(data);
        } else {
            // Another retry needed: cycle back with new dt
            auto retryData = std::make_shared<RetrySequenceData>(
                data->meshes, data->t, newDt, data->iteration + 1, false);
            this->addResult(retryData);
        }
    }

    [[nodiscard]] bool isDone() const { return done_; }

private:
    bool done_ = true;  // Start true: cycle not yet entered
};

/// Custom state manager for the retry loop cycle.
/// Overrides canTerminate() to break the cycle when retries are complete
/// (or when the cycle was never entered).
class RetryLoopStateManager
    : public hh::StateManager<1, RetrySequenceData, RetrySequenceData> {
public:
    RetryLoopStateManager(
        std::shared_ptr<RetryLoopState> const &state,
        std::string const &name)
        : hh::StateManager<1, RetrySequenceData, RetrySequenceData>(
              state, name) {}

    [[nodiscard]] bool canTerminate() const override {
        this->state()->lock();
        auto ret = std::dynamic_pointer_cast<RetryLoopState>(
            this->state())->isDone();
        this->state()->unlock();
        return ret;
    }
};

#endif // CHANGE_TIMESTEP_STATE_H
