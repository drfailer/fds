#ifndef DIV_P1_LATE_KERNEL_TASK_H
#define DIV_P1_LATE_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel kernel task for divergence part 1 late.
/// Extracted from predJoinDivExchange barrier to run per-mesh in parallel
/// after the 2N->N fork-join.
class DivP1LateKernelTask
    : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit DivP1LateKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>(
              "DivP1LateKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_divergence_part_1_late_b(data->nm, data->t, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>>
    copy() override {
        return std::make_shared<DivP1LateKernelTask>(this->numberThreads());
    }
};

#endif // DIV_P1_LATE_KERNEL_TASK_H
