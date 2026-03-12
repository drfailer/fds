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

#endif // VELOCITY_CORRECTOR_KERNEL_TASK_H
