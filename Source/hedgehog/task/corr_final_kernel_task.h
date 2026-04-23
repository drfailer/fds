#ifndef CORR_FINAL_KERNEL_TASK_H
#define CORR_FINAL_KERNEL_TASK_H

#include <hedgehog/hedgehog.h>
#include <sstream>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Corrector-final kernel task: VelocityCorrector + VelocityBCEdges + QRAdd + DivPart2.
///
/// Four input types drive four phases of work:
///   1. MeshData<CorrectorPressure> — velocity corrector per-mesh kernels
///      → emits MeshData<PostVelCorr> to barrier
///   2. MeshData<> — velocity BC edges per-mesh kernels (from barrier)
///      → emits MeshData<> to TimestepTask (via subgraph output)
///   3. MeshData<PostDivJoin> — QR addition + WORK1 copy (from fork2 join, non-CC_IBM)
///      → emits MeshData<DivExch> to DivExchangeTask
///   4. MeshData<DivPart2> — divergence part 2 block kernel (from DivExchangeTask)
///      → emits MeshData<CorrectorPressure> to pressure barrier or subgraph output
///
/// Multi-threaded: all phases run per-mesh in parallel.
class CorrFinalKernelTask
    : public hh::AbstractTask<4,
        MeshData<MeshState::CorrectorPressure>,
        MeshData<>,
        MeshData<MeshState::PostDivJoin>,
        MeshData<MeshState::DivPart2>,
        MeshData<MeshState::PostVelCorr>,
        MeshData<>,
        MeshData<MeshState::DivExch>,
        MeshData<MeshState::CorrectorPressure>> {
public:
    explicit CorrFinalKernelTask(size_t numThreads)
        : hh::AbstractTask<4,
              MeshData<MeshState::CorrectorPressure>,
              MeshData<>,
              MeshData<MeshState::PostDivJoin>,
              MeshData<MeshState::DivPart2>,
              MeshData<MeshState::PostVelCorr>,
              MeshData<>,
              MeshData<MeshState::DivExch>,
              MeshData<MeshState::CorrectorPressure>>(
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

    void execute(std::shared_ptr<MeshData<MeshState::PostDivJoin>> data) override {
        fds_divergence_part_1_add_qr_b(data->nm);
        fds_copy_work1_b_to_work1(data->nm);
        this->addResult(retag<MeshState::DivExch>(data));
    }

    void execute(std::shared_ptr<MeshData<MeshState::DivPart2>> data) override {
        int kbar = fds_get_kbar(data->nm);
        fds_divergence_part_2_block_kernel(data->nm, data->dt, 1, kbar);
        this->addResult(retag<MeshState::CorrectorPressure>(data));
    }

    std::shared_ptr<hh::AbstractTask<4,
        MeshData<MeshState::CorrectorPressure>,
        MeshData<>,
        MeshData<MeshState::PostDivJoin>,
        MeshData<MeshState::DivPart2>,
        MeshData<MeshState::PostVelCorr>,
        MeshData<>,
        MeshData<MeshState::DivExch>,
        MeshData<MeshState::CorrectorPressure>>>
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
            << "Phase 3 (QRAddCopy):\\n"
            << "  DIV_P1_ADD_QR\\n"
            << "  COPY_WORK1_B\\n"
            << "Phase 4 (DivPart2):\\n"
            << "  DIV_P2_BLOCK_KERNEL";
        return oss.str();
    }
};

#endif // CORR_FINAL_KERNEL_TASK_H
