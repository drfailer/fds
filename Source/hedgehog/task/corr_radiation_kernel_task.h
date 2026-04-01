#ifndef CORR_RADIATION_KERNEL_TASK_H
#define CORR_RADIATION_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/corr_radiation_data.h"
#include "../fds_fortran_interface.h"

class CorrRadiationKernelTask
    : public hh::AbstractTask<1, CorrRadiationWork, CorrRadiationWork> {
public:
    explicit CorrRadiationKernelTask(size_t numThreads)
        : hh::AbstractTask<1, CorrRadiationWork, CorrRadiationWork>(
              "CorrRadiationKernel", numThreads) {}

    void execute(std::shared_ptr<CorrRadiationWork> work) override {
        work->radQSumPartial = 0.0;
        work->kfst4SumPartial = 0.0;
        fds_compute_radiation_kernel(
            work->nm, work->t, work->radIter,
            &work->radQSumPartial, &work->kfst4SumPartial);
        fds_cccompute_radiation(work->nm, work->t, work->radIter);
        this->addResult(work);
    }

    std::shared_ptr<hh::AbstractTask<1, CorrRadiationWork, CorrRadiationWork>>
    copy() override {
        return std::make_shared<CorrRadiationKernelTask>(
            this->numberThreads());
    }
};

#endif // CORR_RADIATION_KERNEL_TASK_H
