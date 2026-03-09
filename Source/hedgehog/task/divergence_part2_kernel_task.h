#ifndef DIVERGENCE_PART2_KERNEL_TASK_H
#define DIVERGENCE_PART2_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/divergence_part2_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that calls the thread-safe divergence part 2 kernel.
/// Each thread processes one mesh independently.
class DivergencePart2KernelTask
    : public hh::AbstractTask<1, DivergencePart2Work,
                              DivergencePart2Work> {
public:
    explicit DivergencePart2KernelTask(size_t numThreads)
        : hh::AbstractTask<1, DivergencePart2Work,
                           DivergencePart2Work>(
              "DivPart2Kernel", numThreads) {}

    void execute(
        std::shared_ptr<DivergencePart2Work> work) override {
        fds_divergence_part_2_kernel(work->nm, work->dt);
        this->addResult(work);
    }

    std::shared_ptr<hh::AbstractTask<1, DivergencePart2Work,
                                     DivergencePart2Work>>
    copy() override {
        return std::make_shared<DivergencePart2KernelTask>(
            this->numberThreads());
    }
};

#endif // DIVERGENCE_PART2_KERNEL_TASK_H
