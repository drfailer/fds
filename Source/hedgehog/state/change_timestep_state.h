#ifndef CHANGE_TIMESTEP_STATE_H
#define CHANGE_TIMESTEP_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include <iostream>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// State that implements the CHANGE_TIME_STEP_LOOP from main.f90 (lines 637-760).
///
/// Placed after VelPredictorTask in the graph. Collects all mesh tokens, then checks
/// whether any mesh requires a DT reduction (CHANGE_TIME_STEP_INDEX == -1).
///
/// If retry needed: internally re-runs the predictor sequence (DENSITY through
/// VELOCITY_PREDICTOR) with reduced DT, repeating until no retry is needed.
/// This avoids creating a cycle in the Hedgehog graph.
///
/// The retry re-runs: CC_RESTORE_UVW -> DENSITY -> CC_DENSITY -> MESH_EXCHANGE(1)
///   -> VISCOSITY_BC/VELOCITY_FLUX -> HVAC -> INIT_DIV -> WALL_BC/PARTICLE_MOMENTUM
///   -> DIV_PART_1 -> EXCHANGE_DIV_INFO -> DIV_PART_2 -> PRESSURE_ITERATION
///   -> VELOCITY_PREDICTOR -> STOP_CHECK
///
/// Once no retry is needed, emits all tokens forward.
class ChangeTimeStepState : public hh::AbstractState<1, MeshData, MeshData> {
public:
    ChangeTimeStepState(int nmeshes)
        : hh::AbstractState<1, MeshData, MeshData>(),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // Check if any mesh needs a DT reduction (main.f90 lines 753-758)
            int needRetry = 0;
            double newDt = 0.0;
            fds_check_change_time_step(&needRetry, &newDt);

            while (needRetry) {
                // CFL violation: retry with reduced DT and FIRST_PASS=.FALSE.
                fds_set_first_pass(0);

                for (auto &md : collected_) {
                    md->dt = newDt;
                    md->firstPass = false;
                }

                // Re-run the predictor sequence for all meshes with new DT
                // (main.f90 lines 641-732)
                double t = collected_[0]->t;
                double dt = newDt;

                // DENSITY loop (lines 641-645)
                for (auto &md : collected_) {
                    fds_cc_restore_uvw_unlinked(md->nm);  // only on !FIRST_PASS
                    fds_density(t, dt, md->nm);
                }

                // CC_DENSITY (line 649)
                fds_cc_density(t, dt);

                // MESH_EXCHANGE(1) (line 653)
                fds_mesh_exchange(1);

                // COMPUTE_DIVERGENCE_LOOP (lines 678-686)
                // Note: HVAC_BC_IN is skipped on !FIRST_PASS (line 685)
                for (auto &md : collected_) {
                    fds_set_baroclinic_false(md->nm);
                    fds_viscosity_bc(md->nm, 0);
                    fds_velocity_flux(t, dt, md->nm, 0);
                }

                // HVAC solver (lines 690-694) - HVAC_CALC uses FIRST_PASS=.FALSE.
                fds_hvac_calc(t, dt, 0);  // first=false

                // INITIALIZE_DIVERGENCE_INTEGRALS (line 698)
                fds_initialize_divergence_integrals();

                // WALL_BC + PARTICLE_MOMENTUM + DIVERGENCE_PART_1 (lines 700-704)
                for (auto &md : collected_) {
                    fds_wall_bc(t, dt, md->nm);
                    fds_particle_momentum(dt, md->nm);
                    fds_divergence_part_1(t, dt, md->nm);
                }

                // EXCHANGE_DIVERGENCE_INFO (line 708)
                fds_exchange_divergence_info();

                // DIVERGENCE_PART_2 (lines 712-714)
                for (auto &md : collected_) {
                    fds_divergence_part_2(dt, md->nm);
                }

                // PRESSURE_ITERATION_SCHEME (line 718)
                fds_pressure_iteration(t, dt);

                // Initialize CHANGE_TIME_STEP arrays (lines 722-724)
                fds_init_change_time_step(dt);

                // VELOCITY_PREDICTOR (lines 725-727)
                for (auto &md : collected_) {
                    fds_velocity_predictor(t + dt, dt, md->nm);
                }

                // STOP_CHECK(0) (line 732)
                fds_stop_check_zero();

                // Check for instability
                int stopStatus = fds_get_stop_status();
                if (stopStatus != 0) {
                    break;  // Exit retry loop on instability
                }

                // Check again if retry is needed
                fds_check_change_time_step(&needRetry, &newDt);
            }

            // Emit all tokens forward (no more retry needed)
            for (auto &md : collected_) {
                this->addResult(md);
            }
            collected_.clear();
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

#endif // CHANGE_TIMESTEP_STATE_H
