#ifndef PRESSURE_ITERATION_TASKS_H
#define PRESSURE_ITERATION_TASKS_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel pressure solve kernel task.
///
/// Accepts two input types:
///   - MeshData<SolvePhase>: from PostExchangeRouter (single-process mode)
///   - MeshData<Pressure>: from barrier (MPI mode)
///
/// Outputs MeshData<Pressure> with exchangeRound=1 for post-solve exchange.
///
/// Multi-threaded: each clone processes one mesh independently.
class PressureSolveKernelTask
    : public hh::AbstractTask<2,
          MeshData<MeshState::SolvePhase>,
          MeshData<MeshState::Pressure>,
          MeshData<MeshState::Pressure>> {
public:
    explicit PressureSolveKernelTask(size_t kernelThreads, int presFlag)
        : hh::AbstractTask<2,
              MeshData<MeshState::SolvePhase>,
              MeshData<MeshState::Pressure>,
              MeshData<MeshState::Pressure>>(
              "PressureSolveKernel", kernelThreads),
          presFlag_(presFlag) {}

    void execute(std::shared_ptr<MeshData<MeshState::SolvePhase>> spd) override {
        doWork(retag<MeshState::Pressure>(spd));
    }

    void execute(std::shared_ptr<MeshData<MeshState::Pressure>> md) override {
        doWork(md);
    }

    std::shared_ptr<hh::AbstractTask<2,
        MeshData<MeshState::SolvePhase>,
        MeshData<MeshState::Pressure>,
        MeshData<MeshState::Pressure>>> copy() override {
        return std::make_shared<PressureSolveKernelTask>(
            this->numberThreads(), presFlag_);
    }

    // PRES_FLAG values from GLOBAL_CONSTANTS (cons.f90)
    static constexpr int FFT_PRES_FLAG = 0;
    static constexpr int ULMAT_PRES_FLAG = 3;

private:
    void doWork(std::shared_ptr<MeshData<MeshState::Pressure>> md) {
        if (md->iterate_baroclinic || md->pressure_iterations == 1) {
            if (fds_is_cc_ibm()) {
                fds_cc_match_velocity_flux(md->nm);
            } else {
                fds_match_velocity_flux_kernel(md->nm);
            }
        }
        fds_no_flux_kernel(md->nm, md->dt);
        if (fds_is_cc_ibm()) {
            fds_cc_no_flux(md->dt, md->nm, 0); // FORCE_FLG=FALSE
        }
        if (md->pressure_iterations == 1) {
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
        md->exchangeRound = 1;  // post-solve exchange
        this->addResult(md);
    }

    int presFlag_;
};

#endif // PRESSURE_ITERATION_TASKS_H
