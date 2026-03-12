#ifndef CORR_CONDENS_KERNEL_TASK_H
#define CORR_CONDENS_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel kernel task for corrector condensation.
/// Calls CONDENSATION_EVAPORATION_KERNEL per mesh.
class CorrCondensKernelTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit CorrCondensKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "CorrCondensKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_condensation_kernel(data->nm, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<CorrCondensKernelTask>(
            this->numberThreads());
    }
};

#endif // CORR_CONDENS_KERNEL_TASK_H
