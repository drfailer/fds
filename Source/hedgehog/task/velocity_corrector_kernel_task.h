#ifndef VELOCITY_CORRECTOR_KERNEL_TASK_H
#define VELOCITY_CORRECTOR_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/velocity_corrector_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that executes thread-safe velocity corrector kernels.
///
/// This task is the computational core of the velocity corrector sub-graph.
/// It calls thread-safe Fortran kernels that operate directly on MESHES(NM)
/// without using POINT_TO_MESH, enabling parallel execution of multiple
/// meshes concurrently (numThreads = N for N meshes).
///
/// Kernels called:
/// - VELOCITY_CORRECTOR_KERNEL: Updates velocity field (U = U + FVX*DT)
/// - CHECK_DIVERGENCE_KERNEL: Checks divergence constraints
///
/// Thread-safety: Each thread operates on a different mesh (indexed by nm),
/// with no cross-mesh data access or global state modifications in kernels.
class VelocityCorrectorKernelTask
    : public hh::AbstractTask<1, VelocityCorrectorWork, VelocityCorrectorWork> {
public:
    explicit VelocityCorrectorKernelTask(size_t numThreads)
        : hh::AbstractTask<1, VelocityCorrectorWork, VelocityCorrectorWork>(
              "VelocityCorrectorKernel", numThreads) {}

    void execute(std::shared_ptr<VelocityCorrectorWork> work) override {
        // Call thread-safe kernels directly (no POINT_TO_MESH)
        // These wrappers call VELOCITY_CORRECTOR_KERNEL(MESHES(NM), DT)
        // and CHECK_DIVERGENCE_KERNEL(MESHES(NM)) respectively
        fds_velocity_corrector_kernel(work->nm, work->t, work->dt);
        fds_check_divergence_kernel(work->nm);

        // Pass work token downstream (contains original MeshData)
        this->addResult(work);
    }

    std::shared_ptr<hh::AbstractTask<1, VelocityCorrectorWork, VelocityCorrectorWork>>
    copy() override {
        return std::make_shared<VelocityCorrectorKernelTask>(this->numberThreads());
    }
};

#endif // VELOCITY_CORRECTOR_KERNEL_TASK_H
