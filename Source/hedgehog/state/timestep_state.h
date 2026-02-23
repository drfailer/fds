#ifndef TIMESTEP_STATE_H
#define TIMESTEP_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include <iostream>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// State at the end of the corrector phase (end of one full time step).
/// Collects all mesh tokens, performs global output dumps, checks termination.
/// If T < T_END, adjusts DT, increments ICYC for the next step, sets
/// PREDICTOR=TRUE, and re-emits tokens for the next time step.
///
/// ICYC timing: In the original main.f90, ICYC is incremented at the TOP of
/// MAIN_LOOP (before physics), so during time step N the Fortran ICYC == N.
/// Here, main_hh.cpp sets ICYC=1 before step 1's physics. For subsequent
/// steps, we set ICYC=N+1 after step N's dumps, so that step N+1's physics
/// sees the correct value. The dumps for step N use icyc_ which was already
/// set to N at the start of the step.
class TimestepState : public hh::AbstractState<1, MeshData, MeshData> {
public:
    TimestepState(int nmeshes, double tEnd)
        : hh::AbstractState<1, MeshData, MeshData>(),
          nmeshes_(nmeshes), tEnd_(tEnd) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nmeshes_) {
            double t = collected_[0]->t;
            double dt = collected_[0]->dt;

            // CC_IBM end step (corrector side)
            fds_cc_end_step(t, dt, 0);

            // Set DIAGNOSTICS flag based on current ICYC.
            // icyc_ == N was set before this step's physics ran (by main_hh.cpp
            // for step 1, or by the previous cycle's re-emission for step N>1).
            fds_set_diagnostics(icyc_, t, dt);

            // Global output sequence (matches main.f90 lines 975-996)
            fds_exchange_global_outputs(t, dt);
            fds_update_controls(t, dt);

            // Per-mesh dump (must happen after UPDATE_CONTROLS)
            for (auto &md : collected_) {
                fds_dump_mesh_outputs(md->t, md->dt, md->nm);
            }

            fds_dump_global_outputs(t, dt);
            fds_write_strings(t, dt);
            fds_write_diagnostics(t, dt);

            // Stop check
            fds_stop_check(1, t, dt);

            int stopStatus = fds_get_stop_status();

            // Check termination
            if (t >= tEnd_ || stopStatus != 0) {
                done_ = true;
                // Don't emit - graph will terminate
                collected_.clear();
                return;
            }

            // Prepare next time step
            fds_set_predictor(1);  // PREDICTOR=TRUE
            fds_set_first_pass(1); // FIRST_PASS=TRUE for new CHANGE_TIME_STEP_LOOP

            // Adjust DT based on CFL conditions from the velocity predictor.
            // This replaces the logic at the top of MAIN_LOOP in main.f90:
            //   IF (ALL(CHANGE_TIME_STEP_INDEX==1)) DT = MINVAL(DT_NEW)
            //   IF (ANY(CHANGE_TIME_STEP_INDEX==-1)) DT = MINVAL(DT_NEW)
            //   Clip final time step
            double newDt = fds_adjust_dt(t, dt);

            // Increment ICYC for the NEXT time step, matching the original
            // main.f90 where ICYC is set at the top of MAIN_LOOP before physics.
            // This ensures step N+1's physics sees ICYC == N+1.
            icyc_++;
            fds_set_icyc(icyc_);

            // Re-emit tokens for next predictor step with updated DT
            for (auto &md : collected_) {
                md->phase = 0;  // predictor
                md->dt = newDt; // use CFL-adjusted DT
                md->firstPass = true; // new CHANGE_TIME_STEP_LOOP starts with FIRST_PASS=TRUE
                // T remains as-is (already advanced in corrector phase transition)
                this->addResult(md);
            }
            collected_.clear();
        }
    }

    [[nodiscard]] bool isDone() const { return done_; }

private:
    int nmeshes_;
    double tEnd_;
    int icyc_ = 1;  // Matches main_hh.cpp's initial fds_set_icyc(1)
    bool done_ = false;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

/// Custom state manager that implements canTerminate() for the time-stepping cycle.
class TimestepStateManager : public hh::StateManager<1, MeshData, MeshData> {
public:
    explicit TimestepStateManager(std::shared_ptr<TimestepState> const &state)
        : hh::StateManager<1, MeshData, MeshData>(state, "TimestepLoop") {}

    [[nodiscard]] bool canTerminate() const override {
        this->state()->lock();
        auto ret = std::dynamic_pointer_cast<TimestepState>(this->state())->isDone();
        this->state()->unlock();
        return ret;
    }
};

#endif // TIMESTEP_STATE_H
