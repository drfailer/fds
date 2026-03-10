#ifndef CORR_PARTICLE_KERNEL_TASK_H
#define CORR_PARTICLE_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/corr_particle_data.h"
#include "../fds_fortran_interface.h"

/// Parallel kernel task for corrector particle step.
/// Calls PARTICLE_MOMENTUM_KERNEL per mesh.
class CorrParticleKernelTask : public hh::AbstractTask<1, CorrParticleWork, CorrParticleWork> {
public:
    explicit CorrParticleKernelTask(size_t numThreads)
        : hh::AbstractTask<1, CorrParticleWork, CorrParticleWork>("CorrParticleKernel", numThreads) {}

    void execute(std::shared_ptr<CorrParticleWork> work) override {
        fds_particle_momentum_kernel(work->nm, work->dt);
        this->addResult(work);
    }

    std::shared_ptr<hh::AbstractTask<1, CorrParticleWork, CorrParticleWork>> copy() override {
        return std::make_shared<CorrParticleKernelTask>(this->numberThreads());
    }
};

#endif // CORR_PARTICLE_KERNEL_TASK_H
