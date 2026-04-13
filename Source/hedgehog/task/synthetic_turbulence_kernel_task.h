#ifndef SYNTHETIC_TURBULENCE_KERNEL_TASK_H
#define SYNTHETIC_TURBULENCE_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel kernel task for synthetic turbulence.
/// Extracted from meshExch3SynTurb barrier to run per-mesh in parallel.
class SyntheticTurbulenceKernelTask
    : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit SyntheticTurbulenceKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>(
              "SyntheticTurbulenceKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_synthetic_turbulence_if_enabled(data->dt, data->t, data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>>
    copy() override {
        return std::make_shared<SyntheticTurbulenceKernelTask>(this->numberThreads());
    }
};

#endif // SYNTHETIC_TURBULENCE_KERNEL_TASK_H
