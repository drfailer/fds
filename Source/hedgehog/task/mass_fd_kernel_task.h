#ifndef MASS_FD_KERNEL_TASK_H
#define MASS_FD_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel kernel task for MASS_FINITE_DIFFERENCES only.
/// Used when viscosity is block-decomposed separately (predictor path).
class MassFDKernelTask
    : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit MassFDKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>(
              "MassFDKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_mass_finite_differences_kernel(data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>>
    copy() override {
        return std::make_shared<MassFDKernelTask>(this->numberThreads());
    }
};

/// Parallel kernel task for MASS_FINITE_DIFFERENCES + DENSITY.
/// Used when viscosity is block-decomposed separately (corrector path).
class MassFDDensityKernelTask
    : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit MassFDDensityKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>(
              "MassFDDensityKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_mass_finite_differences_kernel(data->nm);
        fds_density_kernel(data->nm, data->t, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>>
    copy() override {
        return std::make_shared<MassFDDensityKernelTask>(this->numberThreads());
    }
};

#endif // MASS_FD_KERNEL_TASK_H
