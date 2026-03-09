#ifndef VELOCITY_PREDICTOR_KERNEL_TASK_H
#define VELOCITY_PREDICTOR_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/velocity_predictor_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that executes thread-safe velocity predictor kernels.
///
/// This task is the computational core of the velocity predictor sub-graph.
/// It calls thread-safe Fortran kernels that operate directly on MESHES(NM)
/// without using POINT_TO_MESH, enabling parallel execution of multiple
/// meshes concurrently (numThreads = N for N meshes).
///
/// Kernels called:
/// - VELOCITY_PREDICTOR_KERNEL: Predicts velocity field (US = U - DT*(FVX + dH/dx))
/// - CHECK_STABILITY_KERNEL: Computes CFL-limited time step DT_NEW(NM)
///
/// Thread-safety: Each thread operates on a different mesh (indexed by nm),
/// with no cross-mesh data access. DT_NEW(NM) and CHANGE_TIME_STEP_INDEX(NM)
/// are indexed writes (safe). Fortran I/O in CHECK_STABILITY_KERNEL (density
/// clipping warnings) is rare and only fires in extreme conditions.
class VelocityPredictorKernelTask
    : public hh::AbstractTask<1, VelocityPredictorWork, VelocityPredictorWork> {
public:
    explicit VelocityPredictorKernelTask(size_t numThreads)
        : hh::AbstractTask<1, VelocityPredictorWork, VelocityPredictorWork>(
              "VelocityPredictorKernel", numThreads) {}

    void execute(std::shared_ptr<VelocityPredictorWork> work) override {
        // Call thread-safe kernel wrappers directly (no POINT_TO_MESH)
        // Pass t+dt as the time argument, matching the original VelPredictorTask
        fds_velocity_predictor_kernel(work->nm, work->t + work->dt, work->dt);

        // Pass work token downstream (contains original MeshData)
        this->addResult(work);
    }

    std::shared_ptr<hh::AbstractTask<1, VelocityPredictorWork, VelocityPredictorWork>>
    copy() override {
        return std::make_shared<VelocityPredictorKernelTask>(this->numberThreads());
    }
};

#endif // VELOCITY_PREDICTOR_KERNEL_TASK_H
