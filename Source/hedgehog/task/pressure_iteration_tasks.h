#ifndef PRESSURE_ITERATION_TASKS_H
#define PRESSURE_ITERATION_TASKS_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/pressure_iteration_data.h"
#include "../fds_fortran_interface.h"

/// Pre-kernel task for the pressure iteration sub-graph.
///
/// Accepts BarrierData (first entry from outside) or PressureIterData (cycle
/// back from convergence check). Runs pressure iteration init (first entry
/// only), increments the iteration counter, then executes Phase 1:
///   - Baroclinic correction (if needed)
///   - MESH_EXCHANGE(5)
///   - Match velocity flux
///
/// Scatters individual MeshData tokens for parallel Phase 2 kernel.
class PressurePreKernelTask
    : public hh::AbstractTask<2, BarrierData, PressureIterData, MeshData> {
public:
    PressurePreKernelTask()
        : hh::AbstractTask<2, BarrierData, PressureIterData, MeshData>(
              "PressurePreKernel", 1) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        fds_pressure_iteration_init();
        runPhase1(data->meshes, data->t(), data->dt());
    }

    void execute(std::shared_ptr<PressureIterData> data) override {
        runPhase1(data->meshes, data->t, data->dt);
    }

private:
    void runPhase1(std::vector<std::shared_ptr<MeshData>>& meshes,
                   double t, double /*dt*/) {
        fds_pressure_iteration_increment();

        if (fds_pressure_iteration_needs_baroclinic()) {
            for (auto& md : meshes) {
                fds_baroclinic_correction(t, md->nm);
            }
            fds_mesh_exchange(5);
            for (auto& md : meshes) {
                fds_match_velocity_flux_kernel(md->nm);
            }
        }

        // Scatter MeshData for parallel kernel
        for (auto& md : meshes) {
            this->addResult(md);
        }
    }
};

/// Parallel pressure solve kernel task.
///
/// Runs Phase 2 per mesh: NO_FLUX -> WALL_WORK1 zeroing (iteration 1 only) ->
/// COMPUTE_RHS -> solver (FFT or ULMAT) -> CHECK_RESIDUALS.
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
