#ifndef VELOCITY_ERROR_TASK_H
#define VELOCITY_ERROR_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../data/pressure_iteration_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that computes velocity error for a mesh.
///
/// Receives VelErrorPhaseData from PostExchangeRouter (post-solve exchange done).
/// Computes per-mesh velocity error before the convergence check barrier.
class VelocityErrorTask
    : public hh::AbstractTask<1, VelErrorPhaseData, MeshData> {
public:
    VelocityErrorTask(size_t numThreads)
        : hh::AbstractTask<1, VelErrorPhaseData, MeshData>(
              "VelocityError", numThreads) {}

    void execute(std::shared_ptr<VelErrorPhaseData> vepd) override {
        auto data = vepd->mesh;
        fds_compute_velocity_error_kernel(data->nm, data->dt);
        if (fds_is_cc_ibm()) {
            fds_cc_compute_velocity_error(data->dt, data->nm);
        }
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, VelErrorPhaseData, MeshData>> copy() override {
        return std::make_shared<VelocityErrorTask>(this->numberThreads());
    }
};

#endif // VELOCITY_ERROR_TASK_H
