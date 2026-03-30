#ifndef VELOCITY_BC_EDGES_TASK_H
#define VELOCITY_BC_EDGES_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that calls MATCH_VELOCITY_KERNEL + VELOCITY_BC_PREPROCESSING +
/// VELOCITY_BC_PROCESS_EDGES_KERNEL + CC_VELOCITY_BC_TS (if CC_IBM) for one mesh.
/// All routines are thread-safe: explicit M% access, no POINT_TO_MESH.
///
/// In corrector mode (runDevices=true), also calls UPDATE_DEVICES_1_TS.
///
/// @param applyToEstimated 1 for predictor (estimated vars), 0 for corrector (actual vars)
/// @param doIBEdges 1 to process immersed boundary edges, 0 to skip
/// @param runDevices true to call UPDATE_DEVICES_1_TS (corrector only)
class VelocityBCEdgesTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    VelocityBCEdgesTask(size_t numThreads, int applyToEstimated,
                        int doIBEdges = 1, bool runDevices = false)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "VelocityBCEdges", numThreads),
          applyToEstimated_(applyToEstimated), doIBEdges_(doIBEdges),
          runDevices_(runDevices) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_match_velocity_kernel(data->nm, applyToEstimated_);
        fds_velocity_bc_preprocessing(
            data->nm, data->t, applyToEstimated_);
        fds_velocity_bc_process_edges_kernel(
            data->nm, data->t, applyToEstimated_);
        fds_cc_velocity_bc_ts(data->t, data->nm, applyToEstimated_, doIBEdges_);
        if (runDevices_) {
            fds_update_devices_1_ts(data->t, data->dt, data->nm);
        }
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<VelocityBCEdgesTask>(
            this->numberThreads(), applyToEstimated_, doIBEdges_, runDevices_);
    }

private:
    int applyToEstimated_;
    int doIBEdges_;
    bool runDevices_;
};

#endif // VELOCITY_BC_EDGES_TASK_H
