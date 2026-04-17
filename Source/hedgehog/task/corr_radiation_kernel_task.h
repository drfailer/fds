#ifndef CORR_RADIATION_KERNEL_TASK_H
#define CORR_RADIATION_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/corr_radiation_data.h"
#include "../fds_fortran_interface.h"

template<MeshState InS = MeshState::Default>
class CorrRadiationKernelTask
    : public hh::AbstractTask<1, MeshData<InS>, CorrRadiationWork> {
public:
    explicit CorrRadiationKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData<InS>, CorrRadiationWork>(
              "CorrRadiationKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData<InS>> tagged) override {
        auto data = retag<MeshState::Default>(tagged);
        auto work = std::make_shared<CorrRadiationWork>(data->nm, data->t, 1, data);
        fds_compute_radiation_kernel(
            work->nm, work->t, work->radIter,
            &work->radQSumPartial, &work->kfst4SumPartial);
        fds_cccompute_radiation(work->nm, work->t, work->radIter);
        this->addResult(work);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<InS>, CorrRadiationWork>>
    copy() override {
        return std::make_shared<CorrRadiationKernelTask<InS>>(
            this->numberThreads());
    }
};

#endif // CORR_RADIATION_KERNEL_TASK_H
