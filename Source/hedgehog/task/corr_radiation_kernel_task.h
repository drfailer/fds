#ifndef CORR_RADIATION_KERNEL_TASK_H
#define CORR_RADIATION_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include <mutex>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

template<MeshState InS = MeshState::Default>
class CorrRadiationKernelTask
    : public hh::AbstractTask<1, MeshData<InS>, MeshData<MeshState::MeshExch2>> {
public:
    explicit CorrRadiationKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData<InS>, MeshData<MeshState::MeshExch2>>(
              "CorrRadiationKernel", numThreads),
          sumMutex_(std::make_shared<std::mutex>()) {}

    void execute(std::shared_ptr<MeshData<InS>> tagged) override {
        auto data = retag<MeshState::Default>(tagged);
        double radQPartial = 0.0, kfst4Partial = 0.0;
        fds_compute_radiation_kernel(
            data->nm, data->t, 1,
            &radQPartial, &kfst4Partial);
        fds_cccompute_radiation(data->nm, data->t, 1);
        {
            std::lock_guard<std::mutex> lk(*sumMutex_);
            fds_accumulate_rad_sums(radQPartial, kfst4Partial);
        }
        this->addResult(retag<MeshState::MeshExch2>(data));
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<InS>, MeshData<MeshState::MeshExch2>>>
    copy() override {
        auto c = std::make_shared<CorrRadiationKernelTask<InS>>(
            this->numberThreads());
        c->sumMutex_ = this->sumMutex_;
        return c;
    }

private:
    std::shared_ptr<std::mutex> sumMutex_;
};

#endif // CORR_RADIATION_KERNEL_TASK_H
