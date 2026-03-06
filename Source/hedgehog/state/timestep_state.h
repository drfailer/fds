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
/// done flag: if true, the simulation is finished - emits BarrierData to graph
/// output for termination.  If false, re-emits the individual MeshData tokens
/// with updated phase/DT for the next predictor step.
///
/// ICYC timing: TimestepTask increments ICYC and stores the new value in
/// BarrierData::newIcyc.  This state only routes tokens.
class TimestepLoopState : public hh::AbstractState<1, BarrierData, MeshData, BarrierData> {
public:
    TimestepLoopState()
        : hh::AbstractState<1, BarrierData, MeshData, BarrierData>() {}

    void execute(std::shared_ptr<BarrierData> data) override {
        if (data->done) {
            done_ = true;
            // Emit BarrierData to graph output (different type than MeshData cycle)
            this->addResult(data);
            return;
        }

        // Re-emit MeshData tokens for next predictor step with updated DT (cycles back)
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

/// State manager for the time-stepping cycle.
/// No custom canTerminate() needed - the graph terminates when final BarrierData
/// is emitted to the output and all queues drain naturally.
using TimestepLoopStateManager = hh::StateManager<1, BarrierData, MeshData, BarrierData>;

/// Simple sink state for graph termination.
/// Receives the final BarrierData (with done=true) and passes it to graph output.
/// This allows the graph to have a dedicated termination output that doesn't
/// interfere with the MeshData cycle.
class TerminationSinkState : public hh::AbstractState<1, BarrierData, BarrierData> {
public:
    TerminationSinkState()
        : hh::AbstractState<1, BarrierData, BarrierData>() {}

    void execute(std::shared_ptr<BarrierData> data) override {
        this->addResult(data);
    }
};

#endif // TIMESTEP_STATE_H
