#ifndef TIMESTEP_STATE_H
#define TIMESTEP_STATE_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../fds_fortran_interface.h"

/// Pure data-flow state for the time-stepping cycle.
///
/// Receives a BarrierData from TimestepDumpCollector (which has already
/// performed all dump I/O, diagnostics, stop check, and DT adjustment).
/// Checks the done flag: if true, the simulation is finished - emits
/// BarrierData to graph output for termination.  If false, re-emits the
/// individual MeshData tokens with updated phase/DT for the next predictor
/// step.
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

/// Custom state manager for the time-stepping cycle.
/// Overrides canTerminate() to break the cycle when the simulation is done.
/// Without this, Hedgehog cannot terminate nodes in the cycle because they
/// wait on each other indefinitely.
class TimestepLoopStateManager
    : public hh::StateManager<1, BarrierData, MeshData, BarrierData> {
public:
    TimestepLoopStateManager(
        std::shared_ptr<TimestepLoopState> const &state,
        std::string const &name)
        : hh::StateManager<1, BarrierData, MeshData, BarrierData>(
              state, name) {}

    [[nodiscard]] bool canTerminate() const override {
        this->state()->lock();
        auto ret = std::dynamic_pointer_cast<TimestepLoopState>(
            this->state())->isDone();
        this->state()->unlock();
        return ret;
    }
};

/// Collects N MeshData tokens after per-mesh dump I/O, then runs
/// global finalization (DUMP_GLOBAL_OUTPUTS, WRITE_STRINGS, WRITE_DIAGNOSTICS,
/// STOP_CHECK) and decides whether the simulation is done.
///
/// Emits BarrierData with done/newDt/newIcyc for the downstream
/// TimestepLoopState to route.
class TimestepDumpCollector : public hh::AbstractState<1, MeshData, BarrierData> {
public:
    explicit TimestepDumpCollector(int nmeshes, double tEnd,
                                   std::shared_ptr<int> icyc)
        : nmeshes_(nmeshes), tEnd_(tEnd), icyc_(std::move(icyc)),
          nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_[data->nm - nmOffset_] = data;
        ++count_;

        if (count_ == nmeshes_) {
            double t = collected_[0]->t;
            double dt = collected_[0]->dt;

            // Global finalization (must run after all per-mesh dumps)
            fds_dump_global_outputs(t, dt);
            fds_write_strings(t, dt);
            fds_write_diagnostics(t, dt);
            fds_stop_check(1, t, dt);

            auto bd = std::make_shared<BarrierData>();
            int stopStatus = fds_get_stop_status();

            if (t >= tEnd_ || stopStatus != 0) {
                bd->done = true;
            } else {
                bd->done = false;
                fds_set_predictor(1);
                fds_set_first_pass(1);
                bd->newDt = fds_adjust_dt(t, dt);
                (*icyc_)++;
                fds_set_icyc(*icyc_);
                bd->newIcyc = *icyc_;
            }

            bd->meshes = std::move(collected_);
            collected_.resize(nmeshes_, nullptr);
            count_ = 0;
            this->addResult(bd);
        }
    }

private:
    int nmeshes_;
    double tEnd_;
    std::shared_ptr<int> icyc_;
    int nmOffset_;
    int count_ = 0;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

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
