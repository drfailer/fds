#ifndef PARTICLE_MASS_ENERGY_KERNEL_TASK_H
#define PARTICLE_MASS_ENERGY_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel kernel task for particle mass/energy transfer.
/// Calls PARTICLE_MASS_ENERGY_KERNEL per mesh (thread-safe, no POINT_TO_MESH).
class ParticleMassEnergyKernelTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit ParticleMassEnergyKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "ParticleMassEnergy", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_particle_mass_energy_kernel(data->nm, data->t, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<ParticleMassEnergyKernelTask>(
            this->numberThreads());
    }
};

#endif // PARTICLE_MASS_ENERGY_KERNEL_TASK_H
