#ifndef VELOCITY_CORRECTOR_KERNEL_TASK_H
#define VELOCITY_CORRECTOR_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that executes thread-safe velocity corrector kernels.
///
/// Calls per-mesh Fortran kernels matching the full VELOCITY_CORRECTOR
/// subroutine sequence:
///   1. CC_PROJECT_VELOCITY (store, pre-kernel)
///   2. WALL_VELOCITY_NO_GRADH (store, pre-kernel)
///   3. VELOCITY_CORRECTOR_KERNEL
///   4. CC_PROJECT_VELOCITY (fix, post-kernel)
///   5. WALL_VELOCITY_NO_GRADH (fix, post-kernel)
///   6. CHECK_DIVERGENCE_KERNEL
class VelocityCorrectorKernelTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit VelocityCorrectorKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "VelocityCorrectorKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        // Pre-kernel: store wall velocities
        fds_cc_project_velocity_kernel(data->nm, data->dt, 1, 0);
        fds_wall_velocity_no_gradh_kernel(data->nm, data->dt, 1, 0);
        // Main kernel
        fds_velocity_corrector_kernel(data->nm, data->t, data->dt);
        // Post-kernel: fix wall velocities
        fds_cc_project_velocity_kernel(data->nm, data->dt, 0, 0);
        fds_wall_velocity_no_gradh_kernel(data->nm, data->dt, 0, 0);
        // Diagnostic
        fds_check_divergence_kernel(data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<VelocityCorrectorKernelTask>(this->numberThreads());
    }
};

// ============================================================================
// UNUSED — Previously used for non-block fallback paths
// ============================================================================

/// Single-threaded task that calls the full VELOCITY_CORRECTOR subroutine.
/// Must be single-threaded because VELOCITY_CORRECTOR uses POINT_TO_MESH.
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
