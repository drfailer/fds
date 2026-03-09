#ifndef DIV_SETUP_KERNEL_TASK_H
#define DIV_SETUP_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/div_setup_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that calls the thread-safe velocity flux kernel.
/// Each thread processes one mesh independently.
class DivSetupKernelTask
    : public hh::AbstractTask<1, DivSetupWork, DivSetupWork> {
public:
    explicit DivSetupKernelTask(size_t numThreads)
        : hh::AbstractTask<1, DivSetupWork, DivSetupWork>(
              "DivSetupKernel", numThreads) {}

    void execute(std::shared_ptr<DivSetupWork> work) override {
        fds_velocity_flux_kernel(work->nm, work->t, work->dt,
                                 work->estimated);
        this->addResult(work);
    }

    std::shared_ptr<
        hh::AbstractTask<1, DivSetupWork, DivSetupWork>>
    copy() override {
        return std::make_shared<DivSetupKernelTask>(
            this->numberThreads());
    }
};

#endif // DIV_SETUP_KERNEL_TASK_H
