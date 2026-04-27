#ifndef FORKABLE_PARALLEL_COMPUTE_LANE_H
#define FORKABLE_PARALLEL_COMPUTE_LANE_H

#include <hedgehog/hedgehog.h>
#include <atomic>
#include <memory>
#include <mutex>
#include <sstream>
#include <thread_utils/async_worker.hpp>
#include "../data/mesh_data.h"
#include "../data/termination_data.h"
#include "../fds_fortran_interface.h"

/// Unified forkable lane: handles all HEAVY + fork phases across predictor and corrector.
///
/// Phases (input tag → kernels → output tag):
///   P1 (Default):        INSERT_PARTICLES, VISCOSITY, MASS_FD, DENSITY, CC_DENSITY → MeshExch1
///   P2 (PostPredExch):   DIV_P1_PREFORK, Fork[WallBC+DivEarlyB ‖ DivSetup+PartMom], DIV_P1_LATE → DivExch
///   C1 (CorrInput):      VISCOSITY, MASS_FD, DENSITY, CC_DENSITY → MeshExch4
///   C2 (PostCorrStep1):  Fork[Combustion+Soot ‖ DivSetup], emit Default(HVAC), ParticleOps → MeshExch7
///   C3 (PostHvac):       WallBC → PreDivP1 or PostWallBC (HT3D), Radiation → MeshExch2
///
/// Each copy() creates its own AsyncWorker for fork phases P2 and C2.
/// Real OS threads = 2 × numThreads during fork execution.
class ForkableParallelComputeLane
    : public hh::AbstractTask<6,
        MeshData<>,                              // P1
        MeshData<MeshState::PostPredExch>,       // P2
        MeshData<MeshState::CorrInput>,          // C1
        MeshData<MeshState::PostCorrStep1>,      // C2
        MeshData<MeshState::PostHvac>,           // C3
        TerminationData,
        MeshData<MeshState::MeshExch1>,          // P1 out
        MeshData<MeshState::DivExch>,            // P2 out
        MeshData<MeshState::MeshExch4>,          // C1 out
        MeshData<>,                              // C2 out (to HVAC)
        MeshData<MeshState::MeshExch7>,          // C2 out
        MeshData<MeshState::PreDivP1>,           // C3 out (!HT3D: direct pipeline)
        MeshData<MeshState::PostWallBC>,         // C3 out (HT3D: to exchange barrier)
        MeshData<MeshState::MeshExch2>> {        // C3 out

    using TaskBase = hh::AbstractTask<6,
        MeshData<>, MeshData<MeshState::PostPredExch>,
        MeshData<MeshState::CorrInput>, MeshData<MeshState::PostCorrStep1>,
        MeshData<MeshState::PostHvac>, TerminationData,
        MeshData<MeshState::MeshExch1>, MeshData<MeshState::DivExch>,
        MeshData<MeshState::MeshExch4>, MeshData<>,
        MeshData<MeshState::MeshExch7>, MeshData<MeshState::PreDivP1>,
        MeshData<MeshState::PostWallBC>, MeshData<MeshState::MeshExch2>>;

    TU_AsyncWorker worker_{};
    std::shared_ptr<std::atomic<bool>> done_;
    bool ht3d_;
    bool useAsyncWorker_;

    // P2 fork work: WallBC + DivP1EarlyB on AsyncWorker
    static void predForkWork(void *rawData, TU_i64) {
        auto *data = static_cast<MeshData<> *>(rawData);
        fds_wall_bc_preprocessing_kernel(
            data->nm, data->t, data->dt_bc, data->call_ht_1d);
        fds_wall_bc_process_cells_kernel(
            data->nm, data->t, data->dt, data->dt_bc, data->call_ht_1d);
        fds_wall_bc_finalize(data->nm, data->t, data->dt_bc, data->call_ht_1d);
        fds_divergence_part_1_early_b(data->nm, data->t, data->dt);
    }

    // C2 fork work: Combustion + Soot on AsyncWorker
    static void corrForkWork(void *rawData, TU_i64) {
        auto *data = static_cast<MeshData<> *>(rawData);
        fds_combustion_kernel(data->nm, data->t, data->dt);
        fds_soot_oxidation_kernel(data->nm, data->dt);
    }

    void predForkWorkInline(MeshData<> *data) {
        fds_wall_bc_preprocessing_kernel(
            data->nm, data->t, data->dt_bc, data->call_ht_1d);
        fds_wall_bc_process_cells_kernel(
            data->nm, data->t, data->dt, data->dt_bc, data->call_ht_1d);
        fds_wall_bc_finalize(data->nm, data->t, data->dt_bc, data->call_ht_1d);
        fds_divergence_part_1_early_b(data->nm, data->t, data->dt);
    }

    void corrForkWorkInline(MeshData<> *data) {
        fds_combustion_kernel(data->nm, data->t, data->dt);
        fds_soot_oxidation_kernel(data->nm, data->dt);
    }

public:
    explicit ForkableParallelComputeLane(size_t numThreads,
                                         bool ht3d = false,
                                         bool useAsyncWorker = true)
        : TaskBase("ForkableParallelComputeLane", numThreads),
          done_(std::make_shared<std::atomic<bool>>(false)),
          ht3d_(ht3d), useAsyncWorker_(useAsyncWorker) {
        if (useAsyncWorker_) tu_aw_init(&worker_);
    }

    ForkableParallelComputeLane(size_t numThreads,
                                std::shared_ptr<std::atomic<bool>> done,
                                bool ht3d, bool useAsyncWorker)
        : TaskBase("ForkableParallelComputeLane", numThreads),
          done_(std::move(done)), ht3d_(ht3d),
          useAsyncWorker_(useAsyncWorker) {
        if (useAsyncWorker_) tu_aw_init(&worker_);
    }

    ~ForkableParallelComputeLane() override {
        if (useAsyncWorker_) tu_aw_fini(&worker_);
    }

    ForkableParallelComputeLane(ForkableParallelComputeLane const &) = delete;
    ForkableParallelComputeLane &operator=(ForkableParallelComputeLane const &) = delete;

    /// P1: PredStep1 kernels → MeshExch1
    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_set_mesh_predictor(data->nm, 1);
        fds_insert_particles(data->t, data->nm);
        fds_compute_viscosity_kernel(data->nm, 0);
        fds_mass_finite_differences_kernel(data->nm);
        fds_density_kernel(data->nm, data->t, data->dt);
        fds_cc_density_ts(data->nm, data->t, data->dt);
        this->addResult(retag<MeshState::MeshExch1>(data));
    }

    /// P2: Prefork + Fork + DivP1Late → DivExch
    void execute(std::shared_ptr<MeshData<MeshState::PostPredExch>> tagged) override {
        auto data = retag<MeshState::Default>(tagged);
        fds_divergence_part_1_prefork(data->nm, data->t, data->dt);

        if (useAsyncWorker_) {
            tu_aw_exec(&worker_, predForkWork, data.get(), data->nm);
        }

        fds_set_baroclinic_false(data->nm);
        fds_viscosity_bc_kernel(data->nm, data->phase);
        fds_cc_velocity_bc_ts(data->t, data->nm, data->phase, 0);
        fds_velocity_flux_kernel(data->nm, data->t, data->dt, data->phase);
        fds_particle_momentum_kernel(data->nm, data->dt);

        if (useAsyncWorker_) {
            tu_aw_wait(&worker_);
        } else {
            predForkWorkInline(data.get());
        }

        fds_divergence_part_1_late_b(data->nm, data->t, data->dt);
        this->addResult(retag<MeshState::DivExch>(data));
    }

    /// C1: CorrStep1 kernels → MeshExch4
    void execute(std::shared_ptr<MeshData<MeshState::CorrInput>> tagged) override {
        auto data = retag<MeshState::Default>(tagged);
        fds_set_mesh_predictor(data->nm, 0);
        fds_compute_viscosity_kernel(data->nm, 1);
        fds_mass_finite_differences_kernel(data->nm);
        fds_density_kernel(data->nm, data->t, data->dt);
        fds_cc_density_ts(data->nm, data->t, data->dt);
        this->addResult(retag<MeshState::MeshExch4>(data));
    }

    /// C2: DivSetup fork + ParticleOps → Default(HVAC) + MeshExch7
    void execute(std::shared_ptr<MeshData<MeshState::PostCorrStep1>> tagged) override {
        auto data = retag<MeshState::Default>(tagged);

        if (useAsyncWorker_) {
            tu_aw_exec(&worker_, corrForkWork, data.get(), data->nm);
        }

        fds_set_baroclinic_false(data->nm);
        fds_viscosity_bc_kernel(data->nm, data->phase);
        fds_cc_velocity_bc_ts(data->t, data->nm, data->phase, 0);
        fds_velocity_flux_kernel(data->nm, data->t, data->dt, data->phase);
        if (data->phase)
            fds_agglomeration(data->dt, data->nm);

        if (useAsyncWorker_) {
            tu_aw_wait(&worker_);
        } else {
            corrForkWorkInline(data.get());
        }

        // Emit to HVAC barrier (before ParticleOps)
        this->addResult(data);

        // ParticleOps
        fds_condensation_kernel(data->nm, data->dt);
        fds_particle_mass_energy_kernel(data->nm, data->t, data->dt);
        fds_remove_particles(data->t, data->nm);
        fds_move_particles(data->t, data->dt, data->nm);
        fds_particle_momentum_kernel(data->nm, data->dt);

        this->addResult(retag<MeshState::MeshExch7>(data));
    }

    /// C3: WallBC + Radiation → PreDivP1/PostWallBC + MeshExch2
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

        // Emit for div path: direct PreDivP1 (!HT3D) or PostWallBC (HT3D → barrier)
        if (ht3d_) {
            this->addResult(retag<MeshState::PostWallBC>(data));
        } else {
            this->addResult(retag<MeshState::PreDivP1>(data));
        }

        // Radiation compute (overlaps with div path when !HT3D)
        double radQPartial = 0.0, kfst4Partial = 0.0;
        fds_compute_radiation_kernel(data->nm, data->t, 1, &radQPartial, &kfst4Partial);
        fds_cccompute_radiation(data->nm, data->t, 1);
        fds_set_rad_slot(data->nm, radQPartial, kfst4Partial);

        this->addResult(retag<MeshState::MeshExch2>(data));
    }

    void execute(std::shared_ptr<TerminationData>) override {
        done_->store(true);
    }

    [[nodiscard]] bool canTerminate() const override { return done_->load(); }

    std::shared_ptr<TaskBase> copy() override {
        return std::make_shared<ForkableParallelComputeLane>(
            this->numberThreads(), done_, ht3d_, useAsyncWorker_);
    }

    [[nodiscard]] std::string extraPrintingInformation() const override {
        std::ostringstream oss;
        oss << "Threads: " << this->numberThreads();
        if (useAsyncWorker_) oss << " (+1 AsyncWorker each)";
        oss << "\\n"
            << "P1 (PredStep1): INSERT_PART, VISC, MASS_FD, DENS\\n"
            << "P2 (Fork+Div): Fork[WallBC+DivEarly || DivSetup+PartMom]\\n"
            << "C1 (CorrStep1): VISC, MASS_FD, DENS\\n"
            << "C2 (Fork+Part): Fork[Comb+Soot || DivSetup] + ParticleOps\\n"
            << "C3 (WallBC+Rad): WallBC → "
            << (ht3d_ ? "PostWallBC" : "PreDivP1")
            << ", Radiation → MeshExch2";
        return oss.str();
    }
};

#endif // FORKABLE_PARALLEL_COMPUTE_LANE_H
