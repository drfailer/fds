#ifndef WALLBC_FINALIZE_KERNEL_TASK_H
#define WALLBC_FINALIZE_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel kernel task for wall BC finalize.
/// Extracted from groupB barrier to run per-mesh in parallel.
class WallBCFinalizeKernelTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit WallBCFinalizeKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "WallBCFinalizeKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_wall_bc_finalize(data->nm, data->t, data->dt_bc, data->call_ht_1d);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<WallBCFinalizeKernelTask>(this->numberThreads());
    }
};

#endif // WALLBC_FINALIZE_KERNEL_TASK_H
