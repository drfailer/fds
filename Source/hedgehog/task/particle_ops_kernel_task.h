#ifndef PARTICLE_OPS_KERNEL_TASK_H
#define PARTICLE_OPS_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel kernel task for all per-mesh corrector particle operations:
/// condensation, mass/energy transfer, remove, move, and momentum.
/// All routines write only to the current mesh (no cross-mesh writes).
/// Particle cross-mesh transfer is buffered into OMESH send buffers
/// for later MESH_EXCHANGE(7).
class ParticleOpsKernelTask
    : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit ParticleOpsKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>(
              "ParticleOpsKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_condensation_kernel(data->nm, data->dt);
        fds_particle_mass_energy_kernel(data->nm, data->t, data->dt);
        fds_remove_particles(data->t, data->nm);
        fds_move_particles(data->t, data->dt, data->nm);
        fds_particle_momentum_kernel(data->nm, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>>
    copy() override {
        return std::make_shared<ParticleOpsKernelTask>(this->numberThreads());
    }
};

#endif // PARTICLE_OPS_KERNEL_TASK_H
