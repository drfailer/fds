#ifndef VELOCITY_BC_EDGES_TASK_H
#define VELOCITY_BC_EDGES_TASK_H

#include <hedgehog/hedgehog.h>
#include <sstream>
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

/// Merged predictor task: VelocityPredictor + SyntheticTurbulence + VelocityBCEdges.
///
/// Phase 1 (MeshData<PressureTag>, from pressure iteration or DivPart2):
///   VELOCITY_PREDICTOR_KERNEL, CC_PROJECT_VELOCITY, WALL_VELOCITY_NO_GRADH,
///   CHECK_STABILITY → emits MeshData<PostVelPred> (to ChangeTimeStepTask)
///
/// Phase 2 (MeshData<PostPredVelExch>, from post-velocity-exchange barrier):
///   SYNTHETIC_TURBULENCE, CC_VELOCITY_CUTFACES, MATCH_VELOCITY,
///   VELOCITY_BC_PREPROCESSING, VELOCITY_BC_PROCESS_EDGES, CC_VELOCITY_BC
///   → emits MeshData<> (to PhaseTransitionTask)
template<MeshState PressureTag = MeshState::Default>
class PredSynTurbVelBCTask
    : public hh::AbstractTask<2,
        MeshData<PressureTag>,                 // Phase 1: from pressure/DivP2
        MeshData<MeshState::PostPredVelExch>,  // Phase 2: from post-vel-exchange
        MeshData<MeshState::PostVelPred>,      // Phase 1 output
        MeshData<>> {                          // Phase 2 output

    using TaskBase = hh::AbstractTask<2,
        MeshData<PressureTag>, MeshData<MeshState::PostPredVelExch>,
        MeshData<MeshState::PostVelPred>, MeshData<>>;

public:
    explicit PredSynTurbVelBCTask(size_t numThreads)
        : TaskBase("PredSynTurbVelBCKernel", numThreads) {}

    /// Phase 1: VelocityPredictor kernels → PostVelPred
    void execute(std::shared_ptr<MeshData<PressureTag>> data) override {
        fds_velocity_predictor_kernel_only(data->nm, data->dt);
        fds_cc_project_velocity_kernel(data->nm, data->dt, 0, 1);
        fds_wall_velocity_no_gradh_kernel(data->nm, data->dt, 0, 1);
        fds_check_stability_kernel_only(data->nm, data->t + data->dt, data->dt);
        this->addResult(retag<MeshState::PostVelPred>(data));
    }

    /// Phase 2: SyntheticTurbulence + VelocityBCEdges → MeshData<>
    void execute(std::shared_ptr<MeshData<MeshState::PostPredVelExch>> dataIn) override {
        auto data = retag<MeshState::Default>(dataIn);
        fds_synthetic_turbulence_if_enabled(data->dt, data->t, data->nm);
        fds_cc_velocity_cutfaces_ts(data->nm, 1);
        fds_match_velocity_kernel(data->nm, 1);
        fds_velocity_bc_preprocessing(data->nm, data->t, 1);
        fds_velocity_bc_process_edges_kernel(data->nm, data->t, 1);
        fds_cc_velocity_bc_ts(data->t, data->nm, 1, 1);
        this->addResult(data);
    }

    std::shared_ptr<TaskBase> copy() override {
        return std::make_shared<PredSynTurbVelBCTask<PressureTag>>(this->numberThreads());
    }

    [[nodiscard]] std::string extraPrintingInformation() const override {
        std::ostringstream oss;
        oss << "Threads: " << this->numberThreads() << "\\n"
            << "Phase 1 (VelocityPredictor):\\n"
            << "  VEL_PRED_KERNEL\\n"
            << "  CC_PROJECT_VELOCITY\\n"
            << "  WALL_VEL_NO_GRADH\\n"
            << "  CHECK_STABILITY\\n"
            << "Phase 2 (SynTurb+VelBC):\\n"
            << "  SYNTHETIC_TURBULENCE\\n"
            << "  CC_VEL_CUTFACES\\n"
            << "  MATCH_VELOCITY\\n"
            << "  VEL_BC_PREPROC\\n"
            << "  VEL_BC_EDGES\\n"
            << "  CC_VEL_BC";
        return oss.str();
    }
};

#endif // VELOCITY_BC_EDGES_TASK_H
