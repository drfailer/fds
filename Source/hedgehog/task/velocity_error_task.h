#ifndef VELOCITY_ERROR_TASK_H
#define VELOCITY_ERROR_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that computes velocity error for a mesh.
///
/// Called after the post-solve exchange to compute per-mesh velocity error
/// before the convergence check barrier.
class VelocityErrorTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    VelocityErrorTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "VelocityError", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_compute_velocity_error_kernel(data->nm, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>> copy() override {
        return std::make_shared<VelocityErrorTask>(this->numberThreads());
    }
};

#endif // VELOCITY_ERROR_TASK_H
