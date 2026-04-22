#ifndef CORR_DIVSETUP_COMB_PART_TASK_H
#define CORR_DIVSETUP_COMB_PART_TASK_H

#include <hedgehog/hedgehog.h>
#include <sstream>
#include <thread_utils/async_worker.hpp>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Merged corrector task: CorrStep1 + Fork1(DivSetup || Combustion+Soot) + ParticleOps + WallBC.
///
/// Phase 1 (MeshData<>, from subgraph input):
///   CorrStep1 kernels: VISCOSITY, MASS_FD, DENSITY, CC_DENSITY
///   → emits MeshData<MeshExch4> to exchange graph
///
/// Phase 2 (MeshData<PostCorrStep1>, from exchange graph):
///   1. Fork via AsyncWorker:
///      - Worker thread: COMBUSTION + SOOT_OXIDATION
///      - Main thread:   DivSetup (BAROCLINIC, VISC_BC, CC_VEL_BC, VEL_FLUX, AGGLOM)
///   2. Join: tu_aw_wait
///   3. addResult(MeshData<>) → feeds HVAC barrier
///   4. ParticleOps (CONDENSATION, MASS_ENERGY, REMOVE, MOVE, MOMENTUM)
///   5. addResult(MeshData<PostParticleOps>) → feeds MeshExch7 barrier
///
/// Phase 3 (MeshData<PostHvac>, from HvacCalc barrier):
///   WallBC: orch_per_mesh + preprocessing + process_cells + finalize
///   → emits MeshData<PostWallBC> to Fork2 (radiation || DivP1)
///
/// Each copy() creates its own AsyncWorker — no sharing between threads.
/// Real OS thread count = 2 x numThreads (main + worker per thread).
class CorrDivSetupCombPartTask
    : public hh::AbstractTask<3,
        MeshData<>,                              // Phase 1: from subgraph input
        MeshData<MeshState::PostCorrStep1>,      // Phase 2: from exchange graph
        MeshData<MeshState::PostHvac>,           // Phase 3: from HvacCalc barrier
        MeshData<MeshState::MeshExch4>,          // → exchange graph (Phase 1 output)
        MeshData<>,                              // → HVAC barrier (pre-ParticleOps)
        MeshData<MeshState::MeshExch7>,          // → exchange graph (post-ParticleOps)
        MeshData<MeshState::PostWallBC>> {       // → Fork2 (Phase 3 output)

    using TaskBase = hh::AbstractTask<3,
        MeshData<>, MeshData<MeshState::PostCorrStep1>, MeshData<MeshState::PostHvac>,
        MeshData<MeshState::MeshExch4>, MeshData<>,
        MeshData<MeshState::MeshExch7>, MeshData<MeshState::PostWallBC>>;

    TU_AsyncWorker worker_{};

    static void fork1CombWork(void *rawData, TU_i64) {
        auto *data = static_cast<MeshData<> *>(rawData);
        fds_combustion_kernel(data->nm, data->t, data->dt);
        fds_soot_oxidation_kernel(data->nm, data->dt);
    }

public:
    explicit CorrDivSetupCombPartTask(size_t numThreads)
        : TaskBase("CorrDivSetupCombPart", numThreads) {
        tu_aw_init(&worker_);
    }

    ~CorrDivSetupCombPartTask() override { tu_aw_fini(&worker_); }

    CorrDivSetupCombPartTask(CorrDivSetupCombPartTask const &) = delete;
    CorrDivSetupCombPartTask &operator=(CorrDivSetupCombPartTask const &) = delete;

    /// Phase 1: CorrStep1 kernels → emit to exchange graph
    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_compute_viscosity_kernel(data->nm, 1);
        fds_mass_finite_differences_kernel(data->nm);
        fds_density_kernel(data->nm, data->t, data->dt);
        fds_cc_density_ts(data->nm, data->t, data->dt);
        this->addResult(retag<MeshState::MeshExch4>(data));
    }

    /// Phase 2: DivSetup fork + ParticleOps (from exchange graph)
    void execute(std::shared_ptr<MeshData<MeshState::PostCorrStep1>> tagged) override {
        auto data = retag<MeshState::Default>(tagged);

        tu_aw_exec(&worker_, fork1CombWork, data.get(), data->nm);

        fds_set_baroclinic_false(data->nm);
        fds_viscosity_bc_kernel(data->nm, data->phase);
        fds_cc_velocity_bc_ts(data->t, data->nm, data->phase, 0);
        fds_velocity_flux_kernel(data->nm, data->t, data->dt, data->phase);
        if (data->phase)
            fds_agglomeration(data->dt, data->nm);

        tu_aw_wait(&worker_);

        // Emit to HVAC barrier (before ParticleOps)
        this->addResult(data);

        // ParticleOps (concurrent with HVAC barrier collecting)
        fds_condensation_kernel(data->nm, data->dt);
        fds_particle_mass_energy_kernel(data->nm, data->t, data->dt);
        fds_remove_particles(data->t, data->nm);
        fds_move_particles(data->t, data->dt, data->nm);
        fds_particle_momentum_kernel(data->nm, data->dt);

        // Emit to exchange graph (after ParticleOps)
        this->addResult(retag<MeshState::MeshExch7>(data));
    }

    /// Phase 3: WallBC kernels (from HvacCalc barrier)
    void execute(std::shared_ptr<MeshData<MeshState::PostHvac>> tagged) override {
        auto data = retag<MeshState::Default>(tagged);
        double dt_bc = data->dt_bc;
        int call_ht_1d = data->call_ht_1d;
        if (data->phase == 1) {
            fds_wall_bc_orch_per_mesh(data->nm, data->t, data->wall_counter,
                                      &dt_bc, &call_ht_1d);
            data->dt_bc = dt_bc;
            data->call_ht_1d = call_ht_1d;
        }
        fds_wall_bc_preprocessing_kernel(data->nm, data->t, dt_bc, call_ht_1d);
        fds_wall_bc_process_cells_kernel(data->nm, data->t, data->dt, dt_bc, call_ht_1d);
        fds_wall_bc_finalize(data->nm, data->t, dt_bc, call_ht_1d);
        this->addResult(retag<MeshState::PostWallBC>(data));
    }

    std::shared_ptr<TaskBase> copy() override {
        return std::make_shared<CorrDivSetupCombPartTask>(this->numberThreads());
    }

    [[nodiscard]] std::string extraPrintingInformation() const override {
        std::ostringstream oss;
        oss << "Threads: " << this->numberThreads()
            << " (+1 AsyncWorker each)\\n"
            << "Phase 1 (CorrStep1):\\n"
            << "  VISCOSITY, MASS_FD\\n"
            << "  DENSITY, CC_DENSITY\\n"
            << "Phase 2 (DivSetup+Part):\\n"
            << "  Fork(AsyncWorker):\\n"
            << "    A: SET_BARO, VISC_BC\\n"
            << "       VEL_FLUX, AGGLOM\\n"
            << "    B: COMBUSTION, SOOT_OXID\\n"
            << "  -> addResult (HVAC)\\n"
            << "  CONDENSATION\\n"
            << "  PART_MASS_ENERGY\\n"
            << "  REMOVE/MOVE_PARTICLES\\n"
            << "  PARTICLE_MOMENTUM\\n"
            << "Phase 3 (WallBC):\\n"
            << "  WALL_BC_ORCH_PER_MESH\\n"
            << "  PREPROCESSING\\n"
            << "  PROCESS_CELLS\\n"
            << "  FINALIZE";
        return oss.str();
    }
};

#endif // CORR_DIVSETUP_COMB_PART_TASK_H
