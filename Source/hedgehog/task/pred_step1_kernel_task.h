#ifndef PRED_STEP1_KERNEL_TASK_H
#define PRED_STEP1_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/pred_step1_data.h"
#include "../fds_fortran_interface.h"

/// Parallel kernel task for predictor step 1.
/// Calls COMPUTE_VISCOSITY_KERNEL + MASS_FINITE_DIFFERENCES_NEW_KERNEL per mesh.
class PredStep1KernelTask : public hh::AbstractTask<1, PredStep1Work, PredStep1Work> {
public:
    explicit PredStep1KernelTask(size_t numThreads)
        : hh::AbstractTask<1, PredStep1Work, PredStep1Work>("PredStep1Kernel", numThreads) {}

    void execute(std::shared_ptr<PredStep1Work> work) override {
        fds_compute_viscosity_kernel(work->nm, 0);  // estimated=0 for predictor
        fds_mass_finite_differences_kernel(work->nm);
        this->addResult(work);
    }

    std::shared_ptr<hh::AbstractTask<1, PredStep1Work, PredStep1Work>> copy() override {
        return std::make_shared<PredStep1KernelTask>(this->numberThreads());
    }
};

#endif // PRED_STEP1_KERNEL_TASK_H
