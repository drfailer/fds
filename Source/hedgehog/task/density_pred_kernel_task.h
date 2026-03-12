#ifndef DENSITY_PRED_KERNEL_TASK_H
#define DENSITY_PRED_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that calls the thread-safe density kernel.
/// Each thread processes one mesh independently.
class DensityPredKernelTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit DensityPredKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "DensityPredKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_density_kernel(data->nm, data->t, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<DensityPredKernelTask>(
            this->numberThreads());
    }
};

#endif // DENSITY_PRED_KERNEL_TASK_H
