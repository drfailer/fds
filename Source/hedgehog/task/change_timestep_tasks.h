#ifndef CHANGE_TIMESTEP_TASKS_H
#define CHANGE_TIMESTEP_TASKS_H

#include <hedgehog/hedgehog.h>
#include "../data/change_timestep_data.h"
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../fds_fortran_interface.h"

/// All sequential pre-kernel operations for the CFL retry loop.
///
/// Accepts BarrierData from the subgraph entry (first check) and
/// RetrySequenceData from the cycle (subsequent retries).
///
/// When no retry is needed: emits RetrySequenceData(done=true) which
/// bypasses the kernel via type-based routing to RetryLoopState.
///
/// When retry is needed: runs density, mesh exchange, velocity flux,
/// HVAC, divergence init, and wall BC, then scatters MeshData<> tokens
/// for parallel kernel processing.
class RetryPreKernelTask
    : public hh::AbstractTask<2, BarrierData, RetrySequenceData,
                              MeshData<>, RetrySequenceData> {
public:
    RetryPreKernelTask()
        : hh::AbstractTask<2, BarrierData, RetrySequenceData,
                           MeshData<>, RetrySequenceData>("RetryPreKernel", 1) {}

    /// Entry from subgraph input: check if retry is needed.
    void execute(std::shared_ptr<BarrierData> barrier) override {
        fds_stop_check_zero();

        int needRetry = 0;
        double newDt = 0.0;
        fds_check_change_time_step(&needRetry, &newDt);

        if (!needRetry) {
            auto retryData = std::make_shared<RetrySequenceData>(
                barrier->meshes, barrier->t(), barrier->dt(), -1, true);
            this->addResult(retryData);  // bypass to RetryLoopState
            return;
        }

        auto retryData = std::make_shared<RetrySequenceData>(
            barrier->meshes, barrier->t(), newDt, 0, false);
        doPreKernelWork(retryData);
    }

    /// Entry from cycle: always needs retry (RetryLoopState only cycles when needed).
    void execute(std::shared_ptr<RetrySequenceData> data) override {
        doPreKernelWork(data);
    }

private:
    void doPreKernelWork(std::shared_ptr<RetrySequenceData> &data) {
        fds_set_first_pass(0);

        for (auto &md : data->meshes) {
            md->dt = data->dt;
            md->firstPass = false;
        }

        // Density
        for (auto &md : data->meshes) {
            fds_cc_restore_uvw_unlinked(md->nm);
            fds_density(data->t, data->dt, md->nm);
        }

        // CC_DENSITY + mesh exchange
        fds_cc_density(data->t, data->dt);
        fds_mesh_exchange(1);

        // Velocity flux
        for (auto &md : data->meshes) {
            fds_set_baroclinic_false(md->nm);
            fds_viscosity_bc(md->nm, 0);
            fds_velocity_flux(data->t, data->dt, md->nm, 0);
        }

        // HVAC + divergence integrals
        fds_hvac_calc(data->t, data->dt, 0);
        fds_initialize_divergence_integrals();

        // Wall BC (sequential, uses POINT_TO_MESH)
        for (auto &md : data->meshes) {
            fds_wall_bc(data->t, data->dt, md->nm);
        }

        // Scatter MeshData<> for parallel kernel processing
        this->batchAddResult(data->meshes);
    }
};

/// Parallel kernel for particle momentum + divergence part 1 in retry path.
class RetryMomentumDivKernelTask : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit RetryMomentumDivKernelTask(size_t kernelThreads)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>(
              "RetryMomentumDivKernel", kernelThreads) {}

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>>
    copy() override {
        return std::make_shared<RetryMomentumDivKernelTask>(this->numberThreads());
    }

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_particle_momentum_kernel(data->nm, data->dt);
        fds_divergence_part_1_kernel(data->nm, data->t, data->dt);
        this->addResult(data);
    }
};

#endif // CHANGE_TIMESTEP_TASKS_H
