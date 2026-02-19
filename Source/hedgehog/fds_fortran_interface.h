#ifndef FDS_FORTRAN_INTERFACE_H
#define FDS_FORTRAN_INTERFACE_H

// C++ declarations for Fortran ISO_C_BINDING subroutines
// These map to the BIND(C) wrappers in fds_c_interface.f90

extern "C" {

// Input file setup
void fds_set_input_file(const char *fname, int flen);

// Initialization and finalization
void fds_initialize_all(double *t, double *dt, int *nmeshes);
void fds_finalize_all(double t, double dt);

// State setters
void fds_set_predictor(int flag);
void fds_set_icyc(int val);
void fds_get_t_end(double *t_end);
void fds_get_nmeshes(int *nmeshes);
void fds_zero_q_m_dot();
void fds_adjust_dt(double t, double dt, double *dt_out);

// Predictor phase per-mesh subroutines
void fds_insert_particles(double t, int nm);
void fds_compute_viscosity(int nm, int estimated);
void fds_mass_finite_differences(int nm);
void fds_density(double t, double dt, int nm);
void fds_cc_density(double t, double dt);
void fds_viscosity_bc(int nm, int estimated);
void fds_velocity_flux(double t, double dt, int nm, int estimated);
void fds_wall_bc(double t, double dt, int nm);
void fds_particle_momentum(double dt, int nm);
void fds_divergence_part_1(double t, double dt, int nm);
void fds_divergence_part_2(double dt, int nm);
void fds_init_change_time_step(double dt);
void fds_velocity_predictor(double t_plus_dt, double dt, int nm);
void fds_stop_check_zero();
void fds_velocity_corrector(double t, double dt, int nm);
void fds_match_velocity(int nm);
void fds_velocity_bc(double t, int nm, int estimated);

// Corrector phase per-mesh subroutines
void fds_combustion_bc(int nm);
void fds_combustion(double t, double dt);
void fds_condensation(double dt, int nm);
void fds_particle_mass_energy(double t, double dt, int nm);
void fds_move_particles(double t, double dt, int nm);
void fds_compute_radiation(double t, int nm, int rad_iter);
void fds_agglomeration(double dt, int nm);
void fds_cc_end_step(double t, double dt, int diag);
void fds_check_divergence(int nm);

// Output subroutines (per-mesh)
void fds_update_global_outputs(double t, double dt, int nm);
void fds_dump_mesh_outputs(double t, double dt, int nm);

// Global output subroutines (called once after all meshes complete corrector)
void fds_exchange_global_outputs(double t, double dt);
void fds_update_controls(double t, double dt);
void fds_dump_global_outputs(double t, double dt);
void fds_write_strings(double t, double dt);
void fds_write_diagnostics(double t, double dt);
void fds_set_diagnostics(int icyc, double t, double dt);

// Barrier / exchange subroutines
void fds_mesh_exchange(int code);
void fds_post_receives(int code);
void fds_pressure_iteration(double t, double dt);
void fds_initialize_divergence_integrals();
void fds_exchange_divergence_info();
void fds_create_or_remove_obstructions(double t, double dt);
void fds_global_matrix_reassign(int force);
void fds_rte_source_correction();
void fds_stop_check(int end_code, double t, double dt);
void fds_get_stop_status(int *status);
void fds_synthetic_turbulence(double dt, double t, int nm);
void fds_hvac_calc(double t, double dt, int first);
void fds_set_baroclinic_false(int nm);

} // extern "C"

#endif // FDS_FORTRAN_INTERFACE_H
