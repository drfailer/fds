#ifndef VELOCITY_BC_EDGES_TASK_H
#define VELOCITY_BC_EDGES_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Parallel task that calls CC_VELOCITY_CUTFACES_TS (if CC_IBM) +
/// MATCH_VELOCITY_KERNEL + VELOCITY_BC_PREPROCESSING +
/// VELOCITY_BC_PROCESS_EDGES_KERNEL + CC_VELOCITY_BC_TS (if CC_IBM) for one mesh.
/// All routines are thread-safe: explicit M% access, no POINT_TO_MESH.
///
/// In corrector mode (isCorrFinal=true), also calls UPDATE_DEVICES_1_TS,
/// UPDATE_HRR_TS, UPDATE_MASS_TS, UPDATE_FIRE_SPREAD_OUTPUTS_TS.
/// These write to per-mesh indexed arrays; the sequential reduce happens
/// in CorrFinalDumpTask after all meshes complete.
///
/// @param applyToEstimated 1 for predictor (estimated vars), 0 for corrector (actual vars)
/// @param doIBEdges 1 to process immersed boundary edges, 0 to skip
/// @param isCorrFinal true to run corrector-final per-mesh routines
class VelocityBCEdgesTask
    : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    VelocityBCEdgesTask(size_t numThreads, int applyToEstimated,
                        int doIBEdges = 1, bool isCorrFinal = false)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>(
              "VelocityBCEdges", numThreads),
          applyToEstimated_(applyToEstimated), doIBEdges_(doIBEdges),
          isCorrFinal_(isCorrFinal) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_cc_velocity_cutfaces_ts(data->nm, applyToEstimated_);
        fds_match_velocity_kernel(data->nm, applyToEstimated_);
        fds_velocity_bc_preprocessing(
            data->nm, data->t, applyToEstimated_);
        fds_velocity_bc_process_edges_kernel(
            data->nm, data->t, applyToEstimated_);
        fds_cc_velocity_bc_ts(data->t, data->nm, applyToEstimated_, doIBEdges_);
        if (isCorrFinal_) {
            fds_update_devices_1_ts(data->t, data->dt, data->nm);
            fds_update_hrr_ts(data->dt, data->nm);
            fds_update_mass_ts(data->dt, data->nm);
            fds_update_fire_spread_outputs_ts(data->t, data->dt, data->nm);
        }
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>>
    copy() override {
        return std::make_shared<VelocityBCEdgesTask>(
            this->numberThreads(), applyToEstimated_, doIBEdges_, isCorrFinal_);
    }

private:
    int applyToEstimated_;
    int doIBEdges_;
    bool isCorrFinal_;
};

/// Merged predictor-final task: SyntheticTurbulence + VelocityBCEdges.
/// Eliminates the PredFinal sub-graph and standalone SyntheticTurbulenceKernelTask.
///
/// Calls (per mesh):
///   1. SYNTHETIC_TURBULENCE_IF_ENABLED
///   2. CC_VELOCITY_CUTFACES_TS (CC_IBM)
///   3. MATCH_VELOCITY_KERNEL
///   4. VELOCITY_BC_PREPROCESSING
///   5. VELOCITY_BC_PROCESS_EDGES_KERNEL
///   6. CC_VELOCITY_BC_TS (CC_IBM, DO_IBEDGES=TRUE)
class PredSynTurbVelBCTask
    : public hh::AbstractTask<1, MeshData<MeshState::PostPredVelExch>, MeshData<>> {
public:
    explicit PredSynTurbVelBCTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData<MeshState::PostPredVelExch>, MeshData<>>(
              "PredSynTurbVelBCKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData<MeshState::PostPredVelExch>> dataIn) override {
        auto data = retag<MeshState::Default>(dataIn);
        fds_synthetic_turbulence_if_enabled(data->dt, data->t, data->nm);
        fds_cc_velocity_cutfaces_ts(data->nm, 1);  // applyToEstimated=1
        fds_match_velocity_kernel(data->nm, 1);
        fds_velocity_bc_preprocessing(data->nm, data->t, 1);
        fds_velocity_bc_process_edges_kernel(data->nm, data->t, 1);
        fds_cc_velocity_bc_ts(data->t, data->nm, 1, 1);  // applyToEstimated=1, doIBEdges=1
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<MeshState::PostPredVelExch>, MeshData<>>>
    copy() override {
        return std::make_shared<PredSynTurbVelBCTask>(this->numberThreads());
    }
};

#endif // VELOCITY_BC_EDGES_TASK_H
