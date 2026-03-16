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
///
/// Optional pre/post-exchange operations:
/// - ccDensity: run CC_DENSITY(T,DT) before the exchange (after density loops)
/// - ccEndStep: run CC_END_STEP(T,DT) before the exchange (after velocity pred/corr)
/// - initDiv: run INITIALIZE_DIVERGENCE_INTEGRALS after the exchange
class MeshExchangeTask : public hh::AbstractTask<1, BarrierData, MeshData> {
public:
    explicit MeshExchangeTask(int code, bool ccDensity = false,
                              bool ccEndStep = false, bool initDiv = false)
        : hh::AbstractTask<1, BarrierData, MeshData>(
              "MeshExchange(" + std::to_string(code) + ")", 1),
          code_(code), ccDensity_(ccDensity), ccEndStep_(ccEndStep),
          initDiv_(initDiv) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        if (ccDensity_) { fds_cc_density(data->t(), data->dt()); }
        if (ccEndStep_) { fds_cc_end_step(data->t(), data->dt(), 0); }
        fds_mesh_exchange(code_);
        if (initDiv_) { fds_initialize_divergence_integrals(); }
        for (auto &md : data->meshes) { this->addResult(md); }
    }

private:
    int code_;
    bool ccDensity_;
    bool ccEndStep_;
    bool initDiv_;
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

/// Merged COMBUSTION + HVAC barrier task.
/// Eliminates the intermediate collector between Combustion and HVAC.
class CombustionHvacTask : public hh::AbstractTask<1, BarrierData, MeshData> {
public:
    explicit CombustionHvacTask(int first)
        : hh::AbstractTask<1, BarrierData, MeshData>("Combustion+Hvac", 1),
          first_(first) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        fds_combustion(data->t(), data->dt());
        fds_hvac_calc(data->t(), data->dt(), first_);
        for (auto &md : data->meshes) { this->addResult(md); }
    }

private:
    int first_;
};

/// Sequential SOOT_SURFACE_OXIDATION + HVAC_CALC barrier task.
/// Replaces the SOOT loop and HVAC from CombustionHvacTask after parallel combustion.
class SootHvacTask : public hh::AbstractTask<1, BarrierData, MeshData> {
public:
    explicit SootHvacTask(int first)
        : hh::AbstractTask<1, BarrierData, MeshData>("Soot+Hvac", 1),
          first_(first) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        fds_soot_oxidation_loop(data->dt());
        fds_hvac_calc(data->t(), data->dt(), first_);
        for (auto &md : data->meshes) { this->addResult(md); }
    }

private:
    int first_;
};

/// Sequential REMOVE_PARTICLES + MOVE_PARTICLES barrier task.
/// Runs after parallel ParticleMassEnergyKernelTask, before parallel ParticleMomentumKernelTask.
/// REMOVE_PARTICLES writes to OMESH send buffers (cross-mesh), MOVE_PARTICLES has cross-mesh transfer.
class RemoveMoveParticlesTask : public hh::AbstractTask<1, BarrierData, MeshData> {
public:
    RemoveMoveParticlesTask()
        : hh::AbstractTask<1, BarrierData, MeshData>("RemoveMove", 1) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        for (auto &md : data->meshes) {
            fds_remove_particles(md->t, md->nm);
            fds_move_particles(md->t, md->dt, md->nm);
        }
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

/// Merged HVAC + INITIALIZE_DIVERGENCE_INTEGRALS barrier task.
/// Eliminates the intermediate collector between HVAC and InitDiv.
class HvacInitDivTask : public hh::AbstractTask<1, BarrierData, MeshData> {
public:
    explicit HvacInitDivTask(int first)
        : hh::AbstractTask<1, BarrierData, MeshData>("Hvac+InitDiv", 1),
          first_(first) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        fds_hvac_calc(data->t(), data->dt(), first_);
        fds_initialize_divergence_integrals();
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

#endif // BARRIER_TASKS_H
