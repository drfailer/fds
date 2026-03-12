#ifndef CORR_DIV_PART1_KERNEL_TASK_H
#define CORR_DIV_PART1_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/corr_div_part1_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that calls the thread-safe divergence part 1 kernel.
/// Each thread processes one mesh independently.
class CorrDivPart1KernelTask
    : public hh::AbstractTask<1, CorrDivPart1Work, CorrDivPart1Work> {
public:
    explicit CorrDivPart1KernelTask(size_t numThreads)
        : hh::AbstractTask<1, CorrDivPart1Work, CorrDivPart1Work>(
              "CorrDivPart1Kernel", numThreads) {}

    void execute(std::shared_ptr<CorrDivPart1Work> work) override {
        fds_combustion_bc_kernel(work->nm);
        fds_divergence_part_1_kernel(work->nm, work->t, work->dt);
        this->addResult(work);
    }

    std::shared_ptr<
        hh::AbstractTask<1, CorrDivPart1Work, CorrDivPart1Work>>
    copy() override {
        return std::make_shared<CorrDivPart1KernelTask>(
            this->numberThreads());
    }
};

#endif // CORR_DIV_PART1_KERNEL_TASK_H
