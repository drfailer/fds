#ifndef VELOCITY_PREDICTOR_KERNEL_TASK_H
#define VELOCITY_PREDICTOR_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that executes thread-safe velocity predictor kernels.
///
/// Template parameter InS controls the input MeshData state tag.
/// When InS != Default, the task accepts retagged data from the shared
/// pressure subgraph and converts back to MeshData<> for downstream.
///
/// Calls per-mesh Fortran kernels matching the full VELOCITY_PREDICTOR
/// subroutine sequence:
///   1. VELOCITY_PREDICTOR_KERNEL
///   2. CC_PROJECT_VELOCITY (CC_IBM correction)
///   3. WALL_VELOCITY_NO_GRADH (sparse solver wall fixup)
///   4. CHECK_STABILITY (CFL check)
template<MeshState InS = MeshState::Default>
class VelocityPredictorKernelTask
    : public hh::AbstractTask<1, MeshData<InS>, MeshData<>> {
public:
    explicit VelocityPredictorKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData<InS>, MeshData<>>(
              "VelocityPredictorKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData<InS>> data) override {
        fds_velocity_predictor_kernel_only(data->nm, data->dt);
        fds_cc_project_velocity_kernel(data->nm, data->dt, 0, 1);
        fds_wall_velocity_no_gradh_kernel(data->nm, data->dt, 0, 1);
        fds_check_stability_kernel_only(data->nm, data->t + data->dt, data->dt);
        if constexpr (InS == MeshState::Default) {
            this->addResult(data);
        } else {
            this->addResult(data->template retag<MeshState::Default>());
        }
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<InS>, MeshData<>>>
    copy() override {
        return std::make_shared<VelocityPredictorKernelTask<InS>>(this->numberThreads());
    }
};

// ============================================================================
// UNUSED — Previously used for non-block fallback paths
// ============================================================================

/// Single-threaded task that calls the full VELOCITY_PREDICTOR subroutine.
/// Must be single-threaded because VELOCITY_PREDICTOR uses POINT_TO_MESH.
class VelocityPredictorFullTask
    : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    VelocityPredictorFullTask()
        : hh::AbstractTask<1, MeshData<>, MeshData<>>(
              "VelocityPredictorFull", 1) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_velocity_predictor(data->t + data->dt, data->dt, data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>>
    copy() override {
        return std::make_shared<VelocityPredictorFullTask>();
    }
};

#endif // VELOCITY_PREDICTOR_KERNEL_TASK_H
