#ifndef CORR_PARTICLE_KERNEL_TASK_H
#define CORR_PARTICLE_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel kernel task for corrector particle step.
/// Calls PARTICLE_MOMENTUM_KERNEL per mesh.
class CorrParticleKernelTask
    : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit CorrParticleKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>(
              "CorrParticleKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_particle_momentum_kernel(data->nm, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>>
    copy() override {
        return std::make_shared<CorrParticleKernelTask>(
            this->numberThreads());
    }
};

#endif // CORR_PARTICLE_KERNEL_TASK_H
