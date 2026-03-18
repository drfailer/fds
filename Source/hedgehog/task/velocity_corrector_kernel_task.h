#ifndef VELOCITY_CORRECTOR_KERNEL_TASK_H
#define VELOCITY_CORRECTOR_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that executes thread-safe velocity corrector kernels.
///
/// Calls thread-safe Fortran kernels that operate directly on MESHES(NM)
/// without using POINT_TO_MESH, enabling parallel execution of multiple
/// meshes concurrently.
///
/// Kernels called:
/// - VELOCITY_CORRECTOR_KERNEL: Updates velocity field (U = U + FVX*DT)
/// - CHECK_DIVERGENCE_KERNEL: Checks divergence constraints
class VelocityCorrectorKernelTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit VelocityCorrectorKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "VelocityCorrectorKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_velocity_corrector_kernel(data->nm, data->t, data->dt);
        fds_check_divergence_kernel(data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<VelocityCorrectorKernelTask>(this->numberThreads());
    }
};

/// Single-threaded task that calls the full VELOCITY_CORRECTOR subroutine.
///
/// Unlike VelocityCorrectorKernelTask (which only calls the kernel + CHECK_DIVERGENCE),
/// this calls the complete Fortran subroutine which includes:
///   1. WALL_VELOCITY_NO_GRADH(STORE=TRUE) — store wall velocities before kernel
///   2. VELOCITY_CORRECTOR_KERNEL
///   3. WALL_VELOCITY_NO_GRADH(STORE=FALSE) — fix wall velocities after kernel
/// Followed by CHECK_DIVERGENCE_KERNEL for diagnostic output.
///
/// Must be single-threaded because VELOCITY_CORRECTOR uses POINT_TO_MESH
/// which sets global module pointers.
///
/// Used instead of the block-decomposed sub-graph for sparse pressure solvers.
class VelocityCorrectorFullTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    VelocityCorrectorFullTask()
        : hh::AbstractTask<1, MeshData, MeshData>(
              "VelocityCorrectorFull", 1) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_velocity_corrector(data->t, data->dt, data->nm);
        fds_check_divergence_kernel(data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<VelocityCorrectorFullTask>();
    }
};

#endif // VELOCITY_CORRECTOR_KERNEL_TASK_H
