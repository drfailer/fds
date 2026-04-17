#ifndef CORR_FINAL_KERNEL_TASK_H
#define CORR_FINAL_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include <sstream>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../fds_fortran_interface.h"

/// Merged corrector-final kernel task: VelocityCorrector + VelocityBCEdges + RTESourceCorr.
///
/// Three input types drive three phases of work:
///   1. MeshData<CorrectorPressure> — velocity corrector per-mesh kernels
///      → emits MeshData<PostVelCorr> to CorrFinalOrch
///   2. MeshData<> — velocity BC edges per-mesh kernels (from CorrFinalOrch)
///      → emits MeshData<> to CorrFinalDump
///   3. BarrierData — RTE source correction global call (from CorrFinalOrch)
///      → emits BarrierData to CorrFinalDump
///
/// No canTerminate() needed: CorrFinalOrch (the single-threaded cycle partner)
/// receives TerminationData and terminates first, disconnecting from this task.
/// Multi-threaded: phases 1 and 2 run per-mesh in parallel.
/// Phase 3 runs once per timestep on whichever thread picks it up.
class CorrFinalKernelTask
    : public hh::AbstractTask<3,
        MeshData<MeshState::CorrectorPressure>,
        MeshData<>,
        BarrierData,
        MeshData<MeshState::PostVelCorr>,
        MeshData<>,
        BarrierData> {
public:
    explicit CorrFinalKernelTask(size_t numThreads)
        : hh::AbstractTask<3,
              MeshData<MeshState::CorrectorPressure>,
              MeshData<>, BarrierData,
              MeshData<MeshState::PostVelCorr>,
              MeshData<>, BarrierData>(
              "CorrFinalKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData<MeshState::CorrectorPressure>> data) override {
        fds_cc_project_velocity_kernel(data->nm, data->dt, 1, 0);
        fds_wall_velocity_no_gradh_kernel(data->nm, data->dt, 1, 0);
        fds_velocity_corrector_kernel(data->nm, data->t, data->dt);
        fds_cc_project_velocity_kernel(data->nm, data->dt, 0, 0);
        fds_wall_velocity_no_gradh_kernel(data->nm, data->dt, 0, 0);
        fds_check_divergence_kernel(data->nm);
        this->addResult(retag<MeshState::PostVelCorr>(data));
    }

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_cc_velocity_cutfaces_ts(data->nm, 0);
        fds_match_velocity_kernel(data->nm, 0);
        fds_velocity_bc_preprocessing(data->nm, data->t, 0);
        fds_velocity_bc_process_edges_kernel(data->nm, data->t, 0);
        fds_cc_velocity_bc_ts(data->t, data->nm, 0, 1);
        fds_update_devices_1_ts(data->t, data->dt, data->nm);
        fds_update_hrr_ts(data->dt, data->nm);
        fds_update_mass_ts(data->dt, data->nm);
        fds_update_fire_spread_outputs_ts(data->t, data->dt, data->nm);
        this->addResult(data);
    }

    void execute(std::shared_ptr<BarrierData> data) override {
        fds_rte_source_correction();
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<3,
        MeshData<MeshState::CorrectorPressure>,
        MeshData<>, BarrierData,
        MeshData<MeshState::PostVelCorr>,
        MeshData<>, BarrierData>>
    copy() override {
        return std::make_shared<CorrFinalKernelTask>(this->numberThreads());
    }

    [[nodiscard]] std::string extraPrintingInformation() const override {
        std::ostringstream oss;
        oss << "Phase 1 (VelCorr):\\n"
            << "  CC_PROJECT_VEL(store)\\n"
            << "  WALL_VEL_NO_GRADH(store)\\n"
            << "  VEL_CORRECTOR_KERNEL\\n"
            << "  CC_PROJECT_VEL(fix)\\n"
            << "  WALL_VEL_NO_GRADH(fix)\\n"
            << "  CHECK_DIVERGENCE\\n"
            << "Phase 2 (VelBCEdges):\\n"
            << "  CC_VEL_CUTFACES\\n"
            << "  MATCH_VELOCITY\\n"
            << "  VEL_BC_PREPROC+EDGES\\n"
            << "  CC_VEL_BC_TS\\n"
            << "  UPDATE_DEVICES/HRR/MASS\\n"
            << "Phase 3 (RTE):\\n"
            << "  RTE_SOURCE_CORR";
        return oss.str();
    }
};

#endif // CORR_FINAL_KERNEL_TASK_H
