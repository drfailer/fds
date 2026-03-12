#ifndef DIVERGENCE_PART2_KERNEL_TASK_H
#define DIVERGENCE_PART2_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that calls the thread-safe divergence part 2 kernel.
/// Each thread processes one mesh independently.
class DivergencePart2KernelTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit DivergencePart2KernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "DivPart2Kernel", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_divergence_part_2_kernel(data->nm, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<DivergencePart2KernelTask>(
            this->numberThreads());
    }
};

#endif // DIVERGENCE_PART2_KERNEL_TASK_H
