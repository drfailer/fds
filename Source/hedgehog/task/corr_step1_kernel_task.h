#ifndef CORR_STEP1_KERNEL_TASK_H
#define CORR_STEP1_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that calls the three corrector step 1 kernels:
/// COMPUTE_VISCOSITY_KERNEL, MASS_FINITE_DIFFERENCES_NEW_KERNEL,
/// and DENSITY_KERNEL. Each thread processes one mesh independently.
class CorrStep1KernelTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit CorrStep1KernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "CorrStep1Kernel", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        // estimated=true for corrector phase
        fds_compute_viscosity_kernel(data->nm, 1);
        fds_mass_finite_differences_kernel(data->nm);
        fds_density_kernel(data->nm, data->t, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<CorrStep1KernelTask>(
            this->numberThreads());
    }
};

#endif // CORR_STEP1_KERNEL_TASK_H
