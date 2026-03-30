#ifndef TIMESTEP_STATE_H
#define TIMESTEP_STATE_H

#include <hedgehog/hedgehog.h>
#include <chrono>
#include <iomanip>
#include <memory>
#include <sstream>
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
            md->dt_bc = 0.0;        // reset WallBC state (only set in corrector)
            md->call_ht_1d = 0;     // reset WallBC state (only set in corrector)
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

/// Post-dump collector state: collects N MeshData tokens from parallel
/// DumpMeshOutputsTask, then runs global post-dump operations and termination
/// decision.
///
/// Performs:
///   1. Global post-dump ops (DUMP_GLOBAL_OUTPUTS, WRITE_STRINGS, WRITE_DIAGNOSTICS, STOP_CHECK)
///   2. Termination decision (t >= tEnd or STOP_STATUS)
///   3. DT adjustment for next timestep
/// Emits BarrierData with done/newDt/newIcyc for TimestepLoopState.
class PostDumpState : public hh::AbstractState<1, MeshData, BarrierData> {
public:
    PostDumpState(int nmeshes, double tEnd, std::shared_ptr<int> icyc)
        : nmeshes_(nmeshes), tEnd_(tEnd), icyc_(std::move(icyc)),
          nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_[data->nm - nmOffset_] = data;
        if (++count_ < nmeshes_) return;

        auto t0 = std::chrono::steady_clock::now();

        double t = collected_[0]->t;
        double dt = collected_[0]->dt;

        // Global post-dump finalization
        fds_dump_global_outputs(t, dt);
        fds_write_strings(t, dt);
        fds_write_diagnostics(t, dt);
        fds_stop_check(1, t, dt);

        // Build output with termination decision
        auto bd = std::make_shared<BarrierData>();
        int stopStatus = fds_get_stop_status();

        if (t >= tEnd_ || stopStatus != 0) {
            bd->done = true;
        } else {
            bd->done = false;
            fds_set_predictor(1);
            fds_set_first_pass(1);
            bd->newDt = fds_adjust_dt(t, dt);
            ++(*icyc_);
            fds_set_icyc(*icyc_);
            bd->newIcyc = *icyc_;
        }

        bd->meshes = std::move(collected_);
        collected_.resize(nmeshes_, nullptr);
        count_ = 0;

        auto t1 = std::chrono::steady_clock::now();
        totalTime_ += std::chrono::duration<double>(t1 - t0).count();
        ++invocations_;

        this->addResult(bd);
    }

    [[nodiscard]] std::string info() const {
        std::ostringstream oss;
        oss << "DUMP_GLOBAL_OUTPUTS\\n"
            << "WRITE_STRINGS\\n"
            << "WRITE_DIAGNOSTICS\\n"
            << "STOP_CHECK\\n"
            << std::fixed << std::setprecision(3) << totalTime_ << "s"
            << " / " << invocations_ << " calls";
        if (invocations_ > 0)
            oss << " / avg " << std::setprecision(3)
                << (totalTime_ * 1000.0 / invocations_) << "ms";
        return oss.str();
    }

private:
    int nmeshes_;
    double tEnd_;
    std::shared_ptr<int> icyc_;
    int nmOffset_;
    int count_ = 0;
    std::vector<std::shared_ptr<MeshData>> collected_;
    double totalTime_ = 0.0;
    int invocations_ = 0;
};

/// StateManager wrapping PostDumpState, with extraPrintingInformation().
class PostDumpStateManager
    : public hh::StateManager<1, MeshData, BarrierData> {
public:
    PostDumpStateManager(std::shared_ptr<PostDumpState> const &state,
                         std::string const &name)
        : hh::StateManager<1, MeshData, BarrierData>(state, name) {}

    [[nodiscard]] std::string extraPrintingInformation() const override {
        this->state()->lock();
        auto ret = std::dynamic_pointer_cast<PostDumpState>(
            this->state())->info();
        this->state()->unlock();
        return ret;
    }
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
