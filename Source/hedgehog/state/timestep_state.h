#ifndef TIMESTEP_STATE_H
#define TIMESTEP_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"

/// Pure data-flow state for the time-stepping cycle.
///
/// Receives a BarrierData from TimestepTask (which has already performed all
/// computation: outputs, diagnostics, stop check, DT adjustment).  Checks the
/// done flag: if true, the simulation is finished and no tokens are emitted
/// (graph terminates).  If false, re-emits the individual MeshData tokens with
/// updated phase/DT for the next predictor step.
///
/// ICYC timing: TimestepTask increments ICYC and stores the new value in
/// BarrierData::newIcyc.  This state only routes tokens.
class TimestepLoopState : public hh::AbstractState<1, BarrierData, MeshData> {
public:
    TimestepLoopState()
        : hh::AbstractState<1, BarrierData, MeshData>() {}

    void execute(std::shared_ptr<BarrierData> data) override {
        if (data->done) {
            done_ = true;
            return;
        }

        // Re-emit tokens for next predictor step with updated DT
        for (auto &md : data->meshes) {
            md->phase = 0;          // predictor
            md->dt = data->newDt;   // CFL-adjusted DT
            md->firstPass = true;   // new CHANGE_TIME_STEP_LOOP
            this->addResult(md);
        }
    }

    [[nodiscard]] bool isDone() const { return done_; }

private:
    bool done_ = false;
};

/// Custom state manager that implements canTerminate() for the time-stepping cycle.
class TimestepLoopStateManager : public hh::StateManager<1, BarrierData, MeshData> {
public:
    explicit TimestepLoopStateManager(std::shared_ptr<TimestepLoopState> const &state)
        : hh::StateManager<1, BarrierData, MeshData>(state, "TimestepLoop") {}

    [[nodiscard]] bool canTerminate() const override {
        this->state()->lock();
        auto ret = std::dynamic_pointer_cast<TimestepLoopState>(this->state())->isDone();
        this->state()->unlock();
        return ret;
    }
};

#endif // TIMESTEP_STATE_H
