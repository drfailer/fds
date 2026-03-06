#ifndef CHANGE_TIMESTEP_STATE_H
#define CHANGE_TIMESTEP_STATE_H

#include <hedgehog/hedgehog.h>
#include "../data/change_timestep_data.h"
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// State that manages the retry loop.
/// After each retry sequence completion, checks if another retry is needed.
/// Either cycles back (done=false) or marks as complete (done=true).
class RetryLoopState : public hh::AbstractState<1, RetrySequenceData, RetrySequenceData> {
public:
    RetryLoopState() = default;

    void execute(std::shared_ptr<RetrySequenceData> data) override {
        // If already marked done (bypass from CheckRetryTask), pass through
        if (data->done) {
            this->addResult(data);
            return;
        }

        // Check stop status
        int stopStatus = fds_get_stop_status();
        if (stopStatus != 0) {
            // Stop condition met: mark as done and emit
            data->done = true;
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
            this->addResult(data);
        } else {
            // Another retry needed: cycle back with new dt
            auto retryData = std::make_shared<RetrySequenceData>(
                data->meshes, data->t, newDt, data->iteration + 1, false);
            this->addResult(retryData);
        }
    }
};

#endif // CHANGE_TIMESTEP_STATE_H
