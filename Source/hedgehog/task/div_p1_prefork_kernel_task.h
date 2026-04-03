#ifndef DIV_P1_PREFORK_KERNEL_TASK_H
#define DIV_P1_PREFORK_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel kernel task for divergence part 1 prefork.
/// Extracted from meshExch1DivPrefork barrier to run per-mesh in parallel.
class DivP1PreforkKernelTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit DivP1PreforkKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "DivP1PreforkKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_divergence_part_1_prefork(data->nm, data->t, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<DivP1PreforkKernelTask>(this->numberThreads());
    }
};

#endif // DIV_P1_PREFORK_KERNEL_TASK_H
