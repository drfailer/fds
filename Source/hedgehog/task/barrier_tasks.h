#ifndef BARRIER_TASKS_H
#define BARRIER_TASKS_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../fds_fortran_interface.h"

// ---------------------------------------------------------------------------
// Barrier computation tasks.
//
// Each task receives a BarrierData (all mesh tokens collected by a
// CollectorState), performs the global computation, and emits the individual
// MeshData tokens back into the graph.  All run with numThreads=1 because
// the underlying Fortran routines are global (cross-mesh) operations.
// ---------------------------------------------------------------------------

/// MESH_EXCHANGE barrier task — replaces MeshBarrierState.
class MeshExchangeTask : public hh::AbstractTask<1, BarrierData, MeshData> {
public:
    explicit MeshExchangeTask(int code)
        : hh::AbstractTask<1, BarrierData, MeshData>(
              "MeshExchange(" + std::to_string(code) + ")", 1),
          code_(code) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        fds_mesh_exchange(code_);
        for (auto &md : data->meshes) { this->addResult(md); }
    }

private:
    int code_;
};

/// COMBUSTION_LOAD_BALANCED barrier task — replaces CombustionBarrierState.
class CombustionTask : public hh::AbstractTask<1, BarrierData, MeshData> {
public:
    CombustionTask()
        : hh::AbstractTask<1, BarrierData, MeshData>("Combustion", 1) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        fds_combustion(data->t(), data->dt());
        for (auto &md : data->meshes) { this->addResult(md); }
    }
};

/// HVAC_CALC barrier task — replaces HvacBarrierState.
class HvacTask : public hh::AbstractTask<1, BarrierData, MeshData> {
public:
    explicit HvacTask(int first)
        : hh::AbstractTask<1, BarrierData, MeshData>("HvacCalc", 1),
          first_(first) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        fds_hvac_calc(data->t(), data->dt(), first_);
        for (auto &md : data->meshes) { this->addResult(md); }
    }

private:
    int first_;
};

/// PRESSURE_ITERATION_SCHEME barrier task — replaces PressureBarrierState.
class PressureIterationTask : public hh::AbstractTask<1, BarrierData, MeshData> {
public:
    explicit PressureIterationTask(bool predictor = false)
        : hh::AbstractTask<1, BarrierData, MeshData>("PressureIteration", 1),
          predictor_(predictor) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        fds_pressure_iteration(data->t(), data->dt());
        if (predictor_) {
            fds_init_change_time_step(data->dt());
        }
        for (auto &md : data->meshes) { this->addResult(md); }
    }

private:
    bool predictor_;
};

/// INITIALIZE_DIVERGENCE_INTEGRALS barrier task — replaces InitDivIntegralsBarrier.
class InitDivIntegralsTask : public hh::AbstractTask<1, BarrierData, MeshData> {
public:
    InitDivIntegralsTask()
        : hh::AbstractTask<1, BarrierData, MeshData>("InitDivIntegrals", 1) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        fds_initialize_divergence_integrals();
        for (auto &md : data->meshes) { this->addResult(md); }
    }
};

/// EXCHANGE_DIVERGENCE_INFO + RTE barrier task — replaces DivergenceBarrierState.
class DivergenceExchangeTask : public hh::AbstractTask<1, BarrierData, MeshData> {
public:
    explicit DivergenceExchangeTask(bool corrector = false)
        : hh::AbstractTask<1, BarrierData, MeshData>("DivergenceExchange", 1),
          corrector_(corrector) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        fds_exchange_divergence_info();
        if (corrector_) {
            fds_rte_source_correction();
        }
        fds_global_matrix_reassign(0);
        for (auto &md : data->meshes) { this->addResult(md); }
    }

private:
    bool corrector_;
};

/// Phase transition task — replaces PhaseTransitionState.
/// Sets CORRECTOR=TRUE, advances T, zeros arrays, handles obstructions.
class PhaseTransitionTask : public hh::AbstractTask<1, BarrierData, MeshData> {
public:
    PhaseTransitionTask()
        : hh::AbstractTask<1, BarrierData, MeshData>("PhaseTransition", 1) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        double t = data->t();
        double dt = data->dt();

        fds_set_predictor(0);  // CORRECTOR=TRUE, PREDICTOR=FALSE
        t += dt;
        fds_zero_q_m_dot();
        fds_create_or_remove_obstructions(t, dt);
        fds_cc_end_step(t, dt, 0);

        for (auto &md : data->meshes) {
            md->t = t;
            md->phase = 1;  // corrector
            this->addResult(md);
        }
    }
};

/// CHANGE_TIME_STEP_LOOP task — replaces ChangeTimeStepState.
/// Checks CFL compliance; if retry needed, internally re-runs the predictor
/// sequence with reduced DT until no retry is needed.
class ChangeTimeStepTask : public hh::AbstractTask<1, BarrierData, MeshData> {
public:
    ChangeTimeStepTask()
        : hh::AbstractTask<1, BarrierData, MeshData>("ChangeTimeStep", 1) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        fds_stop_check_zero();

        int needRetry = 0;
        double newDt = 0.0;
        fds_check_change_time_step(&needRetry, &newDt);

        while (needRetry) {
            fds_set_first_pass(0);

            for (auto &md : data->meshes) {
                md->dt = newDt;
                md->firstPass = false;
            }

            double t = data->t();
            double dt = newDt;

            for (auto &md : data->meshes) {
                fds_cc_restore_uvw_unlinked(md->nm);
                fds_density(t, dt, md->nm);
            }

            fds_cc_density(t, dt);
            fds_mesh_exchange(1);

            for (auto &md : data->meshes) {
                fds_set_baroclinic_false(md->nm);
                fds_viscosity_bc(md->nm, 0);
                fds_velocity_flux(t, dt, md->nm, 0);
            }

            fds_hvac_calc(t, dt, 0);
            fds_initialize_divergence_integrals();

            for (auto &md : data->meshes) {
                fds_wall_bc(t, dt, md->nm);
                fds_particle_momentum(dt, md->nm);
                fds_divergence_part_1(t, dt, md->nm);
            }

            fds_exchange_divergence_info();

            for (auto &md : data->meshes) {
                fds_divergence_part_2(dt, md->nm);
            }

            fds_pressure_iteration(t, dt);
            fds_init_change_time_step(dt);

            for (auto &md : data->meshes) {
                fds_velocity_predictor(t + dt, dt, md->nm);
            }

            fds_stop_check_zero();

            int stopStatus = fds_get_stop_status();
            if (stopStatus != 0) { break; }

            fds_check_change_time_step(&needRetry, &newDt);
        }

        for (auto &md : data->meshes) { this->addResult(md); }
    }
};

/// Timestep task — computation part of the old TimestepState.
/// Performs global output, diagnostics, stop check, and DT adjustment.
/// Emits BarrierData with done/newDt/newIcyc fields set for the downstream
/// TimestepLoopState to decide whether to terminate or cycle.
class TimestepTask : public hh::AbstractTask<1, BarrierData, BarrierData> {
public:
    explicit TimestepTask(double tEnd)
        : hh::AbstractTask<1, BarrierData, BarrierData>("TimestepCompute", 1),
          tEnd_(tEnd) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        double t = data->t();
        double dt = data->dt();

        fds_cc_end_step(t, dt, 0);
        fds_set_diagnostics(icyc_, t, dt);
        fds_exchange_global_outputs(t, dt);
        fds_update_controls(t, dt);

        for (auto &md : data->meshes) {
            fds_dump_mesh_outputs(md->t, md->dt, md->nm);
        }

        fds_dump_global_outputs(t, dt);
        fds_write_strings(t, dt);
        fds_write_diagnostics(t, dt);
        fds_stop_check(1, t, dt);

        int stopStatus = fds_get_stop_status();

        if (t >= tEnd_ || stopStatus != 0) {
            data->done = true;
        } else {
            data->done = false;
            fds_set_predictor(1);
            fds_set_first_pass(1);
            data->newDt = fds_adjust_dt(t, dt);
            icyc_++;
            fds_set_icyc(icyc_);
            data->newIcyc = icyc_;
        }

        this->addResult(data);
    }

private:
    double tEnd_;
    int icyc_ = 1;
};

#endif // BARRIER_TASKS_H
