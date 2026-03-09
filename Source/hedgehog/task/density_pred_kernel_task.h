#ifndef DENSITY_PRED_KERNEL_TASK_H
#define DENSITY_PRED_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/density_pred_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that calls the thread-safe density kernel.
/// Each thread processes one mesh independently.
class DensityPredKernelTask
    : public hh::AbstractTask<1, DensityPredWork, DensityPredWork> {
public:
    explicit DensityPredKernelTask(size_t numThreads)
        : hh::AbstractTask<1, DensityPredWork, DensityPredWork>(
              "DensityPredKernel", numThreads) {}

    void execute(std::shared_ptr<DensityPredWork> work) override {
        fds_density_kernel(work->nm, work->t, work->dt);
        this->addResult(work);
    }

    std::shared_ptr<
        hh::AbstractTask<1, DensityPredWork, DensityPredWork>>
    copy() override {
        return std::make_shared<DensityPredKernelTask>(
            this->numberThreads());
    }
};

#endif // DENSITY_PRED_KERNEL_TASK_H
