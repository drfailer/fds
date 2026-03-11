#ifndef WALLBC_KERNEL_TASK_H
#define WALLBC_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/wallbc_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that calls WALL_BC_PROCESS_CELLS_KERNEL.
/// Processes ~90% of wall cells without cross-mesh dependencies.
/// Each thread processes one mesh independently.
class WallBCKernelTask
    : public hh::AbstractTask<1, WallBCWork, WallBCWork> {
public:
    explicit WallBCKernelTask(size_t numThreads)
        : hh::AbstractTask<1, WallBCWork, WallBCWork>(
              "WallBCKernel", numThreads) {}

    void execute(std::shared_ptr<WallBCWork> work) override {
        // Call WALL_BC_PROCESS_CELLS_KERNEL - processes cells without
        // HAS_INTERPOLATED_BC or HAS_BACK_MESH flags
        fds_wall_bc_process_cells_kernel(
            work->nm, work->t, work->dt, work->dt_bc, work->call_ht_1d);
        this->addResult(work);
    }

    std::shared_ptr<hh::AbstractTask<1, WallBCWork, WallBCWork>>
    copy() override {
        return std::make_shared<WallBCKernelTask>(
            this->numberThreads());
    }
};

#endif // WALLBC_KERNEL_TASK_H
