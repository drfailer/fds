#ifndef VELOCITY_ERROR_TASK_H
#define VELOCITY_ERROR_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that computes velocity error for a mesh.
///
/// Accepts two input types:
///   - MeshData<VelErrorPhase>: from PostExchangeRouter (single-process mode)
///   - MeshData<Pressure>: from barrier (MPI mode)
///
/// Outputs MeshData<Pressure> to the convergence barrier.
class VelocityErrorTask
    : public hh::AbstractTask<2,
          MeshData<MeshState::VelErrorPhase>,
          MeshData<MeshState::Pressure>,
          MeshData<MeshState::Pressure>> {
public:
    VelocityErrorTask(size_t numThreads)
        : hh::AbstractTask<2,
              MeshData<MeshState::VelErrorPhase>,
              MeshData<MeshState::Pressure>,
              MeshData<MeshState::Pressure>>(
              "VelocityError", numThreads) {}

    void execute(std::shared_ptr<MeshData<MeshState::VelErrorPhase>> vepd) override {
        doWork(vepd->retag<MeshState::Pressure>());
    }

    void execute(std::shared_ptr<MeshData<MeshState::Pressure>> md) override {
        doWork(md);
    }

    std::shared_ptr<hh::AbstractTask<2,
        MeshData<MeshState::VelErrorPhase>,
        MeshData<MeshState::Pressure>,
        MeshData<MeshState::Pressure>>> copy() override {
        return std::make_shared<VelocityErrorTask>(this->numberThreads());
    }

private:
    void doWork(std::shared_ptr<MeshData<MeshState::Pressure>> data) {
        fds_compute_velocity_error_kernel(data->nm, data->dt);
        if (fds_is_cc_ibm()) {
            fds_cc_compute_velocity_error(data->dt, data->nm);
        }
        this->addResult(data);
    }
};

#endif // VELOCITY_ERROR_TASK_H
