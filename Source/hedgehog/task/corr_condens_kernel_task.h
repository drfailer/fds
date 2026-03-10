#ifndef CORR_CONDENS_KERNEL_TASK_H
#define CORR_CONDENS_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/corr_condens_data.h"
#include "../fds_fortran_interface.h"

/// Parallel kernel task for corrector condensation.
/// Calls CONDENSATION_EVAPORATION_KERNEL per mesh.
class CorrCondensKernelTask : public hh::AbstractTask<1, CorrCondensWork, CorrCondensWork> {
public:
    explicit CorrCondensKernelTask(size_t numThreads)
        : hh::AbstractTask<1, CorrCondensWork, CorrCondensWork>("CorrCondensKernel", numThreads) {}

    void execute(std::shared_ptr<CorrCondensWork> work) override {
        fds_condensation_kernel(work->nm, work->dt);
        this->addResult(work);
    }

    std::shared_ptr<hh::AbstractTask<1, CorrCondensWork, CorrCondensWork>> copy() override {
        return std::make_shared<CorrCondensKernelTask>(this->numberThreads());
    }
};

#endif // CORR_CONDENS_KERNEL_TASK_H
