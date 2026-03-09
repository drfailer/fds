#ifndef CHANGE_TIMESTEP_STATE_H
#define CHANGE_TIMESTEP_STATE_H

#include <hedgehog/hedgehog.h>
#include "../data/change_timestep_data.h"
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// State that manages the retry loop cycle.
///
/// Receives RetrySequenceData after each retry sequence pass.
/// Two output types enable Hedgehog type-based routing:
///   - RetrySequenceData → cycles back to retryDensity for another retry
///   - MeshData → exits the subgraph (retry complete or no retry needed)
///
/// Termination is data-driven: the state records the latest (t, dt) from
/// processed tokens and checks t + dt >= tEnd in canTerminate(). This is
/// monotonic — once the simulation time reaches tEnd, the condition stays true.
class RetryLoopState : public hh::AbstractState<1, RetrySequenceData, RetrySequenceData, MeshData> {
public:
    explicit RetryLoopState(double tEnd) : tEnd_(tEnd) {}

    void execute(std::shared_ptr<RetrySequenceData> data) override {
        // Record latest time values from flowing data
        lastT_ = data->t;
        lastDt_ = data->dt;

        if (data->done) {
            // No retry needed (or retries complete): emit MeshData to exit
            for (auto &md : data->meshes) {
                this->addResult(md);
            }
            return;
        }

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
            // Retries complete: emit MeshData to exit
            for (auto &md : data->meshes) {
                this->addResult(md);
            }
        } else {
            // Another retry needed: cycle back with new dt
            auto retryData = std::make_shared<RetrySequenceData>(
                data->meshes, data->t, newDt, data->iteration + 1, false);
            this->addResult(retryData);
        }
    }

    [[nodiscard]] bool reachedEnd() const {
        return lastT_ + lastDt_ >= tEnd_;
    }

private:
    double tEnd_;
    double lastT_ = 0.0;
    double lastDt_ = 0.0;
};

/// Custom state manager for the retry loop cycle.
/// canTerminate() is data-driven: it checks whether the last processed token's
/// t + dt has reached tEnd. This keeps the cycle alive across all time steps
/// and only allows termination when the simulation time reaches the end.
class RetryLoopStateManager
    : public hh::StateManager<1, RetrySequenceData, RetrySequenceData, MeshData> {
public:
    RetryLoopStateManager(
        std::shared_ptr<RetryLoopState> const &state,
        std::string const &name)
        : hh::StateManager<1, RetrySequenceData, RetrySequenceData, MeshData>(
              state, name) {}

    [[nodiscard]] bool canTerminate() const override {
        this->state()->lock();
        auto ret = std::dynamic_pointer_cast<RetryLoopState>(
            this->state())->reachedEnd();
        this->state()->unlock();
        return ret;
    }
};

#endif // CHANGE_TIMESTEP_STATE_H
