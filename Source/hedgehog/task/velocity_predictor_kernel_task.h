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
/// Kernels called:
/// - VELOCITY_PREDICTOR_KERNEL: Predicts velocity field (US = U - DT*(FVX + dH/dx))
/// - CHECK_STABILITY_KERNEL: Computes CFL-limited time step DT_NEW(NM)
class VelocityPredictorKernelTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit VelocityPredictorKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "VelocityPredictorKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        // Pass t+dt as the time argument, matching the original VelPredictorTask
        fds_velocity_predictor_kernel(data->nm, data->t + data->dt, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<VelocityPredictorKernelTask>(this->numberThreads());
    }
};

#endif // VELOCITY_PREDICTOR_KERNEL_TASK_H
