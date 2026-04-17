#ifndef CORR_DIVSETUP_COMB_PART_TASK_H
#define CORR_DIVSETUP_COMB_PART_TASK_H

#include <hedgehog/hedgehog.h>
#include <sstream>
#include <thread_utils/async_worker.hpp>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Merged corrector task: Fork1(DivSetup || Combustion+Soot) + ParticleOps.
///
/// Replaces DivSetupKernelTask + Fork1CombKernelTask + Fork1JoinTask +
/// ParticleOpsKernelTask (4 graph nodes) by using an AsyncWorker for the fork.
///
/// Per-mesh execution:
///   1. Fork via AsyncWorker:
///      - Worker thread: COMBUSTION + SOOT_OXIDATION
///      - Main thread:   DivSetup (BAROCLINIC, VISC_BC, CC_VEL_BC, VEL_FLUX, AGGLOM)
///   2. Join: tu_aw_wait
///   3. addResult(MeshData<>) → feeds HVAC barrier
///   4. ParticleOps (CONDENSATION, MASS_ENERGY, REMOVE, MOVE, MOMENTUM)
///   5. addResult(MeshData<PostParticleOps>) → feeds MeshExch7 barrier
///
/// Two output types enable HVAC to start collecting while ParticleOps runs.
/// Each copy() creates its own AsyncWorker — no sharing between threads.
/// Real OS thread count = 2 x numThreads (main + worker per thread).
class CorrDivSetupCombPartTask
    : public hh::AbstractTask<1, MeshData<>,
        MeshData<>,                              // → HVAC barrier (pre-ParticleOps)
        MeshData<MeshState::PostParticleOps>> {  // → MeshExch7 barrier (post-ParticleOps)
    TU_AsyncWorker worker_{};

    static void fork1CombWork(void *rawData, TU_i64) {
        auto *data = static_cast<MeshData<> *>(rawData);
        fds_combustion_kernel(data->nm, data->t, data->dt);
        fds_soot_oxidation_kernel(data->nm, data->dt);
    }

public:
    explicit CorrDivSetupCombPartTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData<>,
              MeshData<>, MeshData<MeshState::PostParticleOps>>(
              "CorrDivSetupCombPart", numThreads) {
        tu_aw_init(&worker_);
    }

    ~CorrDivSetupCombPartTask() override { tu_aw_fini(&worker_); }

    CorrDivSetupCombPartTask(CorrDivSetupCombPartTask const &) = delete;
    CorrDivSetupCombPartTask &operator=(CorrDivSetupCombPartTask const &) = delete;

    void execute(std::shared_ptr<MeshData<>> data) override {
        // Fork: worker runs Combustion+Soot, main runs DivSetup
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

        // Emit to MeshExch7 barrier (after ParticleOps)
        this->addResult(retag<MeshState::PostParticleOps>(data));
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>,
        MeshData<>, MeshData<MeshState::PostParticleOps>>>
    copy() override {
        return std::make_shared<CorrDivSetupCombPartTask>(this->numberThreads());
    }

    [[nodiscard]] std::string extraPrintingInformation() const override {
        std::ostringstream oss;
        oss << "Fork(AsyncWorker):\\n"
            << "  A: SET_BARO_FALSE, VISC_BC\\n"
            << "     CC_VEL_BC_TS, VEL_FLUX\\n"
            << "     AGGLOMERATION\\n"
            << "  B: COMBUSTION, SOOT_OXID\\n"
            << "-> addResult (HVAC)\\n"
            << "CONDENSATION\\n"
            << "PART_MASS_ENERGY\\n"
            << "REMOVE/MOVE_PARTICLES\\n"
            << "PARTICLE_MOMENTUM";
        return oss.str();
    }
};

#endif // CORR_DIVSETUP_COMB_PART_TASK_H
