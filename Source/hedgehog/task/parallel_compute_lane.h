#ifndef PARALLEL_COMPUTE_LANE_H
#define PARALLEL_COMPUTE_LANE_H

#include <hedgehog/hedgehog.h>
#include <atomic>
#include <memory>
#include <sstream>
#include "../data/mesh_data.h"
#include "../data/termination_data.h"
#include "../fds_fortran_interface.h"

/// Unified parallel lane: handles all MEDIUM/LIGHT phases across predictor and corrector.
///
/// Phases (input tag → kernels → output tag):
///   P3 (DivPart2, phase=0):     DIV_P2_BLOCK → PressureTag
///   P4 (PressureTag):           VEL_PRED, CC_PROJECT, WALL_VEL, CHECK_STABILITY → PostVelPred
///   P5 (PostPredVelExch):       SYNTURB, VEL_CUTFACES, MATCH_VEL, VEL_BC → PredFinalOutput
///   C4 (PreDivP1):              COMBUSTION_BC, DIV_P1_SKIP_QR → PostDivP1
///   C5 (PostDivJoin):           QR_ADD, COPY_WORK1 → DivExch
///   C6 (DivPart2, phase=1):     DIV_P2_BLOCK → CorrectorPressure
///   C7 (CorrectorPressure):     VEL_CORRECTOR, CC_PROJECT×2, WALL_VEL×2, CHECK_DIV → PostVelCorr
///   C8 (PostCorrFinalBarrier):  VEL_BC, UPDATE_DEVICES/HRR/MASS → Default (subgraph output)
///
/// Template parameter PressureTag:
///   PredictorPressure when pressure subgraph active, PredDivP2Out otherwise.
template<MeshState PressureTag = MeshState::PredDivP2Out>
class ParallelComputeLane
    : public hh::AbstractTask<8,
        MeshData<MeshState::DivPart2>,             // P3/C6
        MeshData<PressureTag>,                     // P4
        MeshData<MeshState::PostPredVelExch>,      // P5
        MeshData<MeshState::PreDivP1>,             // C4
        MeshData<MeshState::PostDivJoin>,          // C5
        MeshData<MeshState::CorrectorPressure>,    // C7
        MeshData<MeshState::PostCorrFinalBarrier>, // C8
        TerminationData,
        MeshData<PressureTag>,                     // P3 out
        MeshData<MeshState::PostVelPred>,          // P4 out
        MeshData<MeshState::PredFinalOutput>,      // P5 out
        MeshData<MeshState::PostDivP1>,            // C4 out
        MeshData<MeshState::DivExch>,              // C5 out
        MeshData<MeshState::CorrectorPressure>,    // C6 out
        MeshData<MeshState::PostVelCorr>,          // C7 out
        MeshData<>> {                              // C8 out

    using TaskBase = hh::AbstractTask<8,
        MeshData<MeshState::DivPart2>, MeshData<PressureTag>,
        MeshData<MeshState::PostPredVelExch>, MeshData<MeshState::PreDivP1>,
        MeshData<MeshState::PostDivJoin>, MeshData<MeshState::CorrectorPressure>,
        MeshData<MeshState::PostCorrFinalBarrier>, TerminationData,
        MeshData<PressureTag>, MeshData<MeshState::PostVelPred>,
        MeshData<MeshState::PredFinalOutput>, MeshData<MeshState::PostDivP1>,
        MeshData<MeshState::DivExch>, MeshData<MeshState::CorrectorPressure>,
        MeshData<MeshState::PostVelCorr>, MeshData<>>;

    std::shared_ptr<std::atomic<bool>> done_;

public:
    explicit ParallelComputeLane(size_t numThreads)
        : TaskBase("ParallelComputeLane", numThreads),
          done_(std::make_shared<std::atomic<bool>>(false)) {}

    ParallelComputeLane(size_t numThreads,
                        std::shared_ptr<std::atomic<bool>> done)
        : TaskBase("ParallelComputeLane", numThreads),
          done_(std::move(done)) {}

    /// P3/C6: DivPart2 → PressureTag (pred) or CorrectorPressure (corr)
    void execute(std::shared_ptr<MeshData<MeshState::DivPart2>> data) override {
        int kbar = fds_get_kbar(data->nm);
        fds_divergence_part_2_block_kernel(data->nm, data->dt, 1, kbar);
        if (data->phase == 0) {
            this->addResult(retag<PressureTag>(data));
        } else {
            this->addResult(retag<MeshState::CorrectorPressure>(data));
        }
    }

    /// P4: VelocityPredictor → PostVelPred
    void execute(std::shared_ptr<MeshData<PressureTag>> data) override {
        fds_velocity_predictor_kernel_only(data->nm, data->dt);
        fds_cc_project_velocity_kernel(data->nm, data->dt, 0, 1);
        fds_wall_velocity_no_gradh_kernel(data->nm, data->dt, 0, 1);
        fds_check_stability_kernel_only(data->nm, data->t + data->dt, data->dt);
        this->addResult(retag<MeshState::PostVelPred>(data));
    }

    /// P5: SyntheticTurbulence + VelocityBCEdges → PredFinalOutput
    void execute(std::shared_ptr<MeshData<MeshState::PostPredVelExch>> dataIn) override {
        auto data = retag<MeshState::Default>(dataIn);
        fds_synthetic_turbulence_if_enabled(data->dt, data->t, data->nm);
        fds_cc_velocity_cutfaces_ts(data->nm, 1);
        fds_match_velocity_kernel(data->nm, 1);
        fds_velocity_bc_preprocessing(data->nm, data->t, 1);
        fds_velocity_bc_process_edges_kernel(data->nm, data->t, 1);
        fds_cc_velocity_bc_ts(data->t, data->nm, 1, 1);
        this->addResult(retag<MeshState::PredFinalOutput>(data));
    }

    /// C4: CombustionBC + DivP1 → PostDivP1
    void execute(std::shared_ptr<MeshData<MeshState::PreDivP1>> data) override {
        fds_combustion_bc_kernel(data->nm);
        fds_divergence_part_1_kernel_skip_qr_b(data->nm, data->t, data->dt);
        this->addResult(retag<MeshState::PostDivP1>(data));
    }

    /// C5: QRAdd + CopyWork1 → DivExch
    void execute(std::shared_ptr<MeshData<MeshState::PostDivJoin>> data) override {
        fds_divergence_part_1_add_qr_b(data->nm);
        fds_copy_work1_b_to_work1(data->nm);
        this->addResult(retag<MeshState::DivExch>(data));
    }

    /// C7: VelocityCorrector → PostVelCorr
    void execute(std::shared_ptr<MeshData<MeshState::CorrectorPressure>> data) override {
        fds_cc_project_velocity_kernel(data->nm, data->dt, 1, 0);
        fds_wall_velocity_no_gradh_kernel(data->nm, data->dt, 1, 0);
        fds_velocity_corrector_kernel(data->nm, data->t, data->dt);
        fds_cc_project_velocity_kernel(data->nm, data->dt, 0, 0);
        fds_wall_velocity_no_gradh_kernel(data->nm, data->dt, 0, 0);
        fds_check_divergence_kernel(data->nm);
        this->addResult(retag<MeshState::PostVelCorr>(data));
    }

    /// C8: VelocityBCEdges + UpdateDevices → Default (subgraph output)
    void execute(std::shared_ptr<MeshData<MeshState::PostCorrFinalBarrier>> data) override {
        fds_cc_velocity_cutfaces_ts(data->nm, 0);
        fds_match_velocity_kernel(data->nm, 0);
        fds_velocity_bc_preprocessing(data->nm, data->t, 0);
        fds_velocity_bc_process_edges_kernel(data->nm, data->t, 0);
        fds_cc_velocity_bc_ts(data->t, data->nm, 0, 1);
        fds_update_devices_1_ts(data->t, data->dt, data->nm);
        fds_update_hrr_ts(data->dt, data->nm);
        fds_update_mass_ts(data->dt, data->nm);
        fds_update_fire_spread_outputs_ts(data->t, data->dt, data->nm);
        this->addResult(retag<MeshState::Default>(data));
    }

    void execute(std::shared_ptr<TerminationData>) override {
        done_->store(true);
    }

    [[nodiscard]] bool canTerminate() const override { return done_->load(); }

    std::shared_ptr<TaskBase> copy() override {
        return std::make_shared<ParallelComputeLane<PressureTag>>(
            this->numberThreads(), done_);
    }

    [[nodiscard]] std::string extraPrintingInformation() const override {
        std::ostringstream oss;
        oss << "Threads: " << this->numberThreads() << "\\n"
            << "P3 (DivPart2): DIV_P2_BLOCK\\n"
            << "P4 (VelPred): VEL_PRED, CC_PROJECT, CHECK_STABILITY\\n"
            << "P5 (VelBC): SYNTURB, MATCH_VEL, VEL_BC\\n"
            << "C4 (DivP1): COMBUSTION_BC, DIV_P1_SKIP_QR\\n"
            << "C5 (QRAdd): DIV_P1_ADD_QR, COPY_WORK1\\n"
            << "C6 (DivPart2): DIV_P2_BLOCK\\n"
            << "C7 (VelCorr): VEL_CORRECTOR, CHECK_DIV\\n"
            << "C8 (VelBCEdges): VEL_BC, UPDATE_DEVICES/HRR/MASS";
        return oss.str();
    }
};

#endif // PARALLEL_COMPUTE_LANE_H
