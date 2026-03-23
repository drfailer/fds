#ifndef PARTICLE_MOMENTUM_KERNEL_TASK_H
#define PARTICLE_MOMENTUM_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel kernel task for particle momentum transfer.
/// Calls PARTICLE_MOMENTUM_TRANSFER_KERNEL per mesh (thread-safe).
class ParticleMomentumKernelTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit ParticleMomentumKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "ParticleMomentumKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_particle_momentum_kernel(data->nm, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<ParticleMomentumKernelTask>(
            this->numberThreads());
    }
};

#endif // PARTICLE_MOMENTUM_KERNEL_TASK_H
