#ifndef FDS_FORTRAN_INTERFACE_H
#define FDS_FORTRAN_INTERFACE_H

// C++ declarations for Fortran ISO_C_BINDING subroutines
// These map to the BIND(C) wrappers in fds_c_interface.f90

extern "C" {

// Input file setup
void fds_set_input_file(const char *fname, int flen);

// Mesh re-decomposition: set target mesh dimensions before initialization.
// Each user-configured mesh will be split into sub-meshes of approximately
// these many cells per dimension. Set all to 0 to disable (default).
void fds_set_target_mesh_dims(int ibar, int jbar, int kbar);

// Pressure subgraph override: -1=auto, 0=off, 1=on
void fds_set_pressure_subgraph(int mode);

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
int fds_get_kbar(int nm);
int fds_is_cc_ibm();
int fds_has_radiation();
int fds_exchange_radiation();
int fds_use_pressure_subgraph();
int fds_iterate_pressure();
int fds_get_pres_flag();
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
void fds_velocity_predictor_block_kernel(int nm, double dt, int k1, int k2);
void fds_check_stability_kernel_only(int nm, double t, double dt);
void fds_velocity_corrector_kernel(int nm, double t, double dt);
void fds_velocity_corrector_block_kernel(int nm, double dt, int k1, int k2);
void fds_check_divergence_kernel(int nm);
void fds_divergence_part_2_kernel(int nm, double dt);
void fds_divergence_part_2_preprocessing(int nm, double dt);
void fds_divergence_part_2_block_kernel(int nm, double dt, int k1, int k2);
int fds_divergence_part_2_can_block_decompose();
void fds_compute_viscosity_kernel(int nm, int estimated);
void fds_compute_viscosity_block_kernel(int nm, int estimated, int k1, int k2);
void fds_compute_viscosity_post_block(int nm, int estimated);
int fds_compute_viscosity_can_block_decompose();
void fds_cutface_velocities(int nm, int estimated, int cutfaces);
void fds_cc_velocity_flux_post(int nm, double t, double dt, int estimated);
void fds_mass_finite_differences_kernel(int nm);
void fds_density_kernel(int nm, double t, double dt);
void fds_density_block_preprocessing(int nm, double t, double dt);
void fds_density_block_kernel(int nm, double t, double dt, int k1, int k2);
void fds_density_block_postprocessing(int nm, double t, double dt);
int fds_density_can_block_decompose();
void fds_divergence_part_1_kernel(int nm, double t, double dt);
void fds_divergence_part_1_kernel_skip_qr(int nm, double t, double dt);
void fds_divergence_part_1_kernel_skip_qr_b(int nm, double t, double dt);
void fds_divergence_part_1_add_qr(int nm);
void fds_divergence_part_1_add_qr_b(int nm);
void fds_divergence_part_1_prefork(int nm, double t, double dt);
void fds_divergence_part_1_early_b(int nm, double t, double dt);
void fds_divergence_part_1_late_b(int nm, double t, double dt);
void fds_velocity_flux_kernel(int nm, double t, double dt, int estimated);
void fds_velocity_flux_block_kernel(int nm, double t, double dt, int estimated, int k1, int k2);
int fds_velocity_flux_can_block_decompose(int nm);
void fds_particle_momentum_kernel(int nm, double dt);
void fds_particle_momentum_block_kernel(int nm, double dt, int k1, int k2);
void fds_viscosity_bc_kernel(int nm, int estimated);
void fds_combustion_bc_kernel(int nm);
void fds_condensation_kernel(int nm, double dt);
void fds_wall_bc_preprocessing(int nm, double t, double dt_bc, int call_ht_1d);
void fds_wall_bc_preprocessing_kernel(int nm, double t, double dt_bc, int call_ht_1d);
void fds_wall_bc_process_cells_kernel(int nm, double t, double dt, double dt_bc, int call_ht_1d);
void fds_wall_bc_process_cells_block_kernel(int nm, double t, double dt, double dt_bc, int call_ht_1d, int k1, int k2);
int fds_wall_bc_can_block_decompose();
void fds_wall_bc_finalize(int nm, double t, double dt_bc, int call_ht_1d);
void fds_match_velocity_kernel(int nm, int is_predictor);
void fds_velocity_bc_preprocessing(int nm, double t, int apply_to_estimated);
void fds_velocity_bc_process_edges_kernel(int nm, double t, int apply_to_estimated);
void fds_velocity_bc_process_edges_block_kernel(int nm, double t, int apply_to_estimated,
                                                 int k1, int k2, double *drag_uvwmax_out);
void fds_set_drag_uvwmax(int nm, double val);
void fds_synthetic_turbulence_if_enabled(double dt, double t, int nm);
void fds_cc_velocity_bc(double t, int nm, int estimated, int do_ibedges);
void fds_cc_project_velocity_kernel(int nm, double dt, int store_flag, int predictor_flag);
void fds_wall_velocity_no_gradh_kernel(int nm, double dt, int store_flag, int predictor_flag);

// Corrector phase per-mesh subroutines
void fds_combustion_bc(int nm);
void fds_combustion(double t, double dt);
void fds_combustion_kernel(int nm, double t, double dt);
void fds_soot_oxidation_loop(double dt);
void fds_condensation(double dt, int nm);
void fds_particle_mass_energy(double t, double dt, int nm);
void fds_particle_mass_energy_kernel(int nm, double t, double dt);
void fds_remove_particles(double t, int nm);
void fds_move_particles(double t, double dt, int nm);
void fds_compute_radiation(double t, int nm, int rad_iter);
void fds_compute_radiation_kernel(int nm, double t, int rad_iter,
    double* rad_q_sum_out, double* kfst4_sum_out);
void fds_compute_radiation_kernel_b(int nm, double t, int rad_iter,
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
void fds_exchange_inserted_particles();
void fds_pressure_iteration(double t, double dt);

// Pressure iteration kernel routines (for sub-graph parallelization)
void fds_no_flux_kernel(int nm, double dt);
void fds_match_velocity_flux_kernel(int nm);
void fds_pressure_solver_compute_rhs_kernel(int nm, double t, double dt);
void fds_pressure_solver_fft_kernel(int nm);
void fds_pressure_check_residuals_kernel(int nm);
void fds_ulmat_solver_kernel(int nm, double t, double dt);
void fds_ulmat_check_residuals_kernel(int nm);
void fds_compute_velocity_error_kernel(int nm, double dt);

// Pressure iteration sub-graph helper functions
void fds_pressure_iteration_init();
void fds_pressure_iteration_increment();
int fds_get_pressure_iterations();
void fds_baroclinic_correction(double t, int nm);
void fds_pressure_iteration_zero_wall_work1(int nm);
void fds_pressure_iteration_check_convergence(double t, double dt);
int fds_pressure_iteration_converged();
int fds_pressure_iteration_needs_baroclinic();

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

// Per-neighbor flux exchange (CODE 5 decomposition)
int fds_flux_get_neighbor_count(int nm);
int fds_flux_get_neighbor_mesh(int nm, int idx);
int fds_flux_has_send_cells(int nm, int nom);
int fds_flux_recv_count(int nm);
void fds_flux_copy_neighbor(int nm, int nom);
void fds_flux_copy_neighbor_ts(int nm, int nom);

// Cross-process flux exchange (pack/unpack for CommunicatorTask)
void fds_flux_pack(int nm, int nom, double *buf, int bufsize);
void fds_flux_unpack(int nm, int nom, const double *buf, int bufsize);
int fds_flux_pack_size(int nm, int nom);
int fds_flux_get_process(int nm);
int fds_flux_max_buffer_size();

// Generic mesh exchange dependency queries
int fds_exchange_recv_dep_count(int nm);   // how many meshes send TO nm
int fds_exchange_recv_dep_mesh(int nm, int idx); // 1-based idx -> NOM
int fds_exchange_send_dep_count(int nm);   // how many meshes nm sends TO
int fds_exchange_send_dep_mesh(int nm, int idx); // 1-based idx -> NOM
int fds_mesh_process(int nm);              // MPI rank owning mesh nm
int fds_get_total_meshes();                // global mesh count
void fds_dump_mesh_exchange_topology();    // debug: dump NIC topology

} // extern "C"

#endif // FDS_FORTRAN_INTERFACE_H
