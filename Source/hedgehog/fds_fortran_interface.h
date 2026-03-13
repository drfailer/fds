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
double fds_get_t_end();
int fds_get_nmeshes();
int fds_get_lower_mesh_index();
int fds_get_upper_mesh_index();
int fds_is_cc_ibm();
void fds_zero_q_m_dot();
double fds_adjust_dt(double t, double dt);

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
void fds_set_first_pass(int flag);
void fds_check_change_time_step(int *need_retry, double *dt_out);
void fds_cc_restore_uvw_unlinked(int nm);
void fds_velocity_predictor(double t_plus_dt, double dt, int nm);
void fds_stop_check_zero();
void fds_velocity_corrector(double t, double dt, int nm);
void fds_match_velocity(int nm);
void fds_velocity_bc(double t, int nm, int estimated);

// Thread-safe kernel wrappers (bypass orchestration, call kernels directly)
void fds_velocity_predictor_kernel(int nm, double t, double dt);
void fds_velocity_predictor_kernel_only(int nm, double dt);
void fds_check_stability_kernel_only(int nm, double t, double dt);
void fds_velocity_corrector_kernel(int nm, double t, double dt);
void fds_check_divergence_kernel(int nm);
void fds_divergence_part_2_kernel(int nm, double dt);
void fds_compute_viscosity_kernel(int nm, int estimated);
void fds_mass_finite_differences_kernel(int nm);
void fds_density_kernel(int nm, double t, double dt);
void fds_divergence_part_1_kernel(int nm, double t, double dt);
void fds_velocity_flux_kernel(int nm, double t, double dt, int estimated);
void fds_particle_momentum_kernel(int nm, double dt);
void fds_viscosity_bc_kernel(int nm, int estimated);
void fds_combustion_bc_kernel(int nm);
void fds_condensation_kernel(int nm, double dt);
void fds_wall_bc_preprocessing(int nm, double t, double dt_bc, int call_ht_1d);
void fds_wall_bc_preprocessing_kernel(int nm, double t, double dt_bc, int call_ht_1d);
void fds_wall_bc_process_cells_kernel(int nm, double t, double dt, double dt_bc, int call_ht_1d);
void fds_wall_bc_finalize(int nm, double t, double dt_bc, int call_ht_1d);
void fds_match_velocity_kernel(int nm, int is_predictor);
void fds_velocity_bc_preprocessing(int nm, double t, int apply_to_estimated);
void fds_velocity_bc_process_edges_kernel(int nm, double t, int apply_to_estimated);
void fds_synthetic_turbulence_if_enabled(double dt, double t, int nm);
void fds_cc_velocity_bc(double t, int nm, int estimated);
void fds_cc_project_velocity(int nm, double dt, int store_flag);
void fds_wall_velocity_no_gradh(int nm, double dt, int store_flag);

// Corrector phase per-mesh subroutines
void fds_combustion_bc(int nm);
void fds_combustion(double t, double dt);
void fds_condensation(double dt, int nm);
void fds_particle_mass_energy(double t, double dt, int nm);
void fds_move_particles(double t, double dt, int nm);
void fds_compute_radiation(double t, int nm, int rad_iter);
void fds_compute_radiation_kernel(int nm, double t, int rad_iter,
    double* rad_q_sum_out, double* kfst4_sum_out);
void fds_accumulate_rad_sums(double rad_q_partial, double kfst4_partial);
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
void fds_flush_output_files();

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
int fds_get_stop_status();
void fds_synthetic_turbulence(double dt, double t, int nm);
void fds_hvac_calc(double t, double dt, int first);
void fds_set_baroclinic_false(int nm);

// WALL_BC helper functions
double fds_compute_wall_bc_dt_bc(double t);
void fds_increment_wall_counter();
int fds_check_call_ht_1d();
void fds_reset_wall_counter();
void fds_update_bc_clock(double t);

} // extern "C"

#endif // FDS_FORTRAN_INTERFACE_H
