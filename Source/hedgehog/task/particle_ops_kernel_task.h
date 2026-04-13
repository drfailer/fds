#ifndef PARTICLE_OPS_KERNEL_TASK_H
#define PARTICLE_OPS_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel kernel task for condensation + particle mass/energy transfer.
/// Extracted from groupA barrier to run per-mesh in parallel.
///
/// Note: fds_remove_particles and fds_move_particles are NOT included
/// because they still use POINT_TO_MESH (not thread-safe). Those remain
/// in the downstream sequential barrier.
class ParticleOpsKernelTask
    : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit ParticleOpsKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>(
              "ParticleOpsKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_condensation_kernel(data->nm, data->dt);
        fds_particle_mass_energy_kernel(data->nm, data->t, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>>
    copy() override {
        return std::make_shared<ParticleOpsKernelTask>(this->numberThreads());
    }
};

#endif // PARTICLE_OPS_KERNEL_TASK_H
