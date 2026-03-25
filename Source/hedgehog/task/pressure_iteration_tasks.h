#ifndef PRESSURE_ITERATION_TASKS_H
#define PRESSURE_ITERATION_TASKS_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel pressure solve kernel task.
///
/// Runs Phase 2 per mesh: MATCH_VELOCITY_FLUX (if baroclinic, moved from
/// PreKernel for parallelization) -> NO_FLUX -> WALL_WORK1 zeroing
/// (iteration 1 only) -> COMPUTE_RHS -> solver (FFT or ULMAT) ->
/// CHECK_RESIDUALS.
///
/// Multi-threaded: each clone processes one mesh independently.
class PressureSolveKernelTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit PressureSolveKernelTask(size_t kernelThreads, int presFlag)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "PressureSolveKernel", kernelThreads),
          presFlag_(presFlag) {}

    void execute(std::shared_ptr<MeshData> md) override {
        // match_velocity_flux_kernel: only called when baroclinic term is active
        // (or first iteration). Guards must match the original Fortran conditional:
        // IF (ITERATE_BAROCLINIC_TERM .OR. PRESSURE_ITERATIONS==1).
        // The flag is stable during parallel execution (set/cleared in barriers).
        if (fds_pressure_iteration_needs_baroclinic() ||
            fds_get_pressure_iterations() == 1) {
            fds_match_velocity_flux_kernel(md->nm);
        }
        fds_no_flux_kernel(md->nm, md->dt);
        if (fds_get_pressure_iterations() == 1) {
            fds_pressure_iteration_zero_wall_work1(md->nm);
        }
        fds_pressure_solver_compute_rhs_kernel(md->nm, md->t, md->dt);
        if (presFlag_ == ULMAT_PRES_FLAG) {
            fds_ulmat_solver_kernel(md->nm, md->t, md->dt);
            fds_ulmat_check_residuals_kernel(md->nm);
        } else {
            fds_pressure_solver_fft_kernel(md->nm);
            fds_pressure_check_residuals_kernel(md->nm);
        }
        this->addResult(md);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>> copy() override {
        return std::make_shared<PressureSolveKernelTask>(
            this->numberThreads(), presFlag_);
    }

    // PRES_FLAG values from GLOBAL_CONSTANTS (cons.f90)
    static constexpr int FFT_PRES_FLAG = 0;
    static constexpr int ULMAT_PRES_FLAG = 3;

private:
    int presFlag_;
};

#endif // PRESSURE_ITERATION_TASKS_H
