#ifndef VELOCITY_BC_EDGES_TASK_H
#define VELOCITY_BC_EDGES_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that calls MATCH_VELOCITY_KERNEL + VELOCITY_BC_PREPROCESSING +
/// VELOCITY_BC_PROCESS_EDGES_KERNEL for one mesh.
/// All three routines are thread-safe: explicit M% access, no POINT_TO_MESH.
///
/// @param applyToEstimated 1 for predictor (estimated vars), 0 for corrector (actual vars)
class VelocityBCEdgesTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    VelocityBCEdgesTask(size_t numThreads, int applyToEstimated)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "VelocityBCEdges", numThreads),
          applyToEstimated_(applyToEstimated) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_match_velocity_kernel(data->nm, applyToEstimated_);
        fds_velocity_bc_preprocessing(
            data->nm, data->t, applyToEstimated_);
        fds_velocity_bc_process_edges_kernel(
            data->nm, data->t, applyToEstimated_);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<VelocityBCEdgesTask>(this->numberThreads(), applyToEstimated_);
    }

private:
    int applyToEstimated_;
};

#endif // VELOCITY_BC_EDGES_TASK_H
