#ifndef VELOCITY_PREDICTOR_KERNEL_TASK_H
#define VELOCITY_PREDICTOR_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that executes thread-safe velocity predictor kernels.
///
/// Calls thread-safe Fortran kernels that operate directly on MESHES(NM)
/// without using POINT_TO_MESH, enabling parallel execution of multiple
/// meshes concurrently.
///
/// For non-CC_IBM cases (skipCFL=false):
/// - VELOCITY_PREDICTOR_KERNEL + CHECK_STABILITY_KERNEL (combined)
///
/// For CC_IBM cases (skipCFL=true):
/// - VELOCITY_PREDICTOR_KERNEL only (CFL check runs later in the CC collector,
///   after CC_PROJECT_VELOCITY and WALL_VELOCITY_NO_GRADH)
class VelocityPredictorKernelTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    VelocityPredictorKernelTask(size_t numThreads, bool skipCFL = false)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "VelocityPredictorKernel", numThreads),
          skipCFL_(skipCFL) {}

    void execute(std::shared_ptr<MeshData> data) override {
        if (skipCFL_) {
            fds_velocity_predictor_kernel_only(data->nm, data->dt);
        } else {
            fds_velocity_predictor_kernel(data->nm, data->t + data->dt, data->dt);
        }
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<VelocityPredictorKernelTask>(this->numberThreads(), skipCFL_);
    }

private:
    bool skipCFL_;
};

/// Single-threaded task that calls the full VELOCITY_PREDICTOR subroutine.
///
/// Unlike VelocityPredictorKernelTask (which only calls the kernel + CFL check),
/// this calls the complete Fortran subroutine which includes:
///   1. VELOCITY_PREDICTOR_KERNEL
///   2. WALL_VELOCITY_NO_GRADH (for ULMAT/GLMAT/UGLMAT sparse pressure solvers)
///   3. CHECK_STABILITY
///
/// Must be single-threaded because VELOCITY_PREDICTOR uses POINT_TO_MESH
/// which sets global module pointers.
///
/// Used instead of the block-decomposed sub-graph for sparse pressure solvers.
class VelocityPredictorFullTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    VelocityPredictorFullTask()
        : hh::AbstractTask<1, MeshData, MeshData>(
              "VelocityPredictorFull", 1) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_velocity_predictor(data->t + data->dt, data->dt, data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<VelocityPredictorFullTask>();
    }
};

#endif // VELOCITY_PREDICTOR_KERNEL_TASK_H
