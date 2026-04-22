#ifndef PRED_PREFORK_DIV_TASK_H
#define PRED_PREFORK_DIV_TASK_H

#include <hedgehog/hedgehog.h>
#include <sstream>
#include <thread_utils/async_worker.hpp>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Merged predictor task: PredStep1 + DivP1Prefork + Fork(DivSetup+PartMom || WallBC+DivP1Early)
/// + DivP1Late + DivPart2.
///
/// Phase 1 (MeshData<>, from subgraph input):
///   INSERT_PARTICLES, VISCOSITY, MASS_FD, DENSITY, CC_DENSITY
///   → emits MeshData<MeshExch1> to exchange graph
///
/// Phase 2 (MeshData<PostPredExch>, from post-exchange barrier):
///   1. DIV_P1_PREFORK (main thread, sequential)
///   2. Fork via AsyncWorker:
///      - Worker thread: WALL_BC + DIV_P1_EARLY_B
///      - Main thread:   DivSetup + PARTICLE_MOMENTUM
///   3. Join: tu_aw_wait
///   4. DIV_P1_LATE_B → emits MeshData<DivExch>
///
/// Phase 3 (MeshData<DivPart2>, from DivExchangeTask):
///   DIV_P2_BLOCK_KERNEL → emits MeshData<PressureTag>
///
/// Each copy() creates its own AsyncWorker — no sharing between Hedgehog
/// threads.  Real OS threads = 2 x numThreads (main + worker per thread).
template<MeshState PressureTag = MeshState::Default>
class PredPreforkDivTask
    : public hh::AbstractTask<3,
        MeshData<>,                        // Phase 1: PredStep1 (from subgraph input)
        MeshData<MeshState::PostPredExch>, // Phase 2: from post-exchange barrier
        MeshData<MeshState::DivPart2>,     // Phase 3: from DivExchangeTask
        MeshData<MeshState::MeshExch1>,    // Phase 1 output (to exchange graph)
        MeshData<MeshState::DivExch>,      // Phase 2 output (to DivExchangeTask)
        MeshData<PressureTag>> {           // Phase 3 output (downstream)

    using TaskBase = hh::AbstractTask<3,
        MeshData<>, MeshData<MeshState::PostPredExch>, MeshData<MeshState::DivPart2>,
        MeshData<MeshState::MeshExch1>, MeshData<MeshState::DivExch>,
        MeshData<PressureTag>>;

    TU_AsyncWorker worker_{};

    static void wallBCDivEarlyWork(void *rawData, TU_i64) {
        auto *data = static_cast<MeshData<> *>(rawData);
        fds_wall_bc_preprocessing_kernel(
            data->nm, data->t, data->dt_bc, data->call_ht_1d);
        fds_wall_bc_process_cells_kernel(
            data->nm, data->t, data->dt, data->dt_bc, data->call_ht_1d);
        fds_wall_bc_finalize(data->nm, data->t, data->dt_bc, data->call_ht_1d);
        fds_divergence_part_1_early_b(data->nm, data->t, data->dt);
    }

public:
    explicit PredPreforkDivTask(size_t numThreads)
        : TaskBase("PredPreforkDiv", numThreads) {
        tu_aw_init(&worker_);
    }

    ~PredPreforkDivTask() override { tu_aw_fini(&worker_); }

    PredPreforkDivTask(PredPreforkDivTask const &) = delete;
    PredPreforkDivTask &operator=(PredPreforkDivTask const &) = delete;

    /// Phase 1: PredStep1 kernels → emit to exchange graph
    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_insert_particles(data->t, data->nm);
        fds_compute_viscosity_kernel(data->nm, 0);
        fds_mass_finite_differences_kernel(data->nm);
        fds_density_kernel(data->nm, data->t, data->dt);
        fds_cc_density_ts(data->nm, data->t, data->dt);
        this->addResult(retag<MeshState::MeshExch1>(data));
    }

    /// Phase 2: Prefork + Fork + DivP1Late (from post-exchange barrier)
    void execute(std::shared_ptr<MeshData<MeshState::PostPredExch>> tagged) override {
        auto data = retag<MeshState::Default>(tagged);
        fds_divergence_part_1_prefork(data->nm, data->t, data->dt);

        tu_aw_exec(&worker_, wallBCDivEarlyWork, data.get(), data->nm);

        fds_set_baroclinic_false(data->nm);
        fds_viscosity_bc_kernel(data->nm, data->phase);
        fds_cc_velocity_bc_ts(data->t, data->nm, data->phase, 0);
        fds_velocity_flux_kernel(data->nm, data->t, data->dt, data->phase);
        fds_particle_momentum_kernel(data->nm, data->dt);

        tu_aw_wait(&worker_);

        fds_divergence_part_1_late_b(data->nm, data->t, data->dt);
        this->addResult(retag<MeshState::DivExch>(data));
    }

    /// Phase 3: DivPart2 block (from DivExchangeTask)
    void execute(std::shared_ptr<MeshData<MeshState::DivPart2>> data) override {
        int kbar = fds_get_kbar(data->nm);
        fds_divergence_part_2_block_kernel(data->nm, data->dt, 1, kbar);
        this->addResult(retag<PressureTag>(data));
    }

    std::shared_ptr<TaskBase> copy() override {
        return std::make_shared<PredPreforkDivTask<PressureTag>>(this->numberThreads());
    }

    [[nodiscard]] std::string extraPrintingInformation() const override {
        std::ostringstream oss;
        oss << "Threads: " << this->numberThreads()
            << " (+1 AsyncWorker each)\\n"
            << "Phase 1 (PredStep1):\\n"
            << "  INSERT_PARTICLES\\n"
            << "  VISCOSITY, MASS_FD\\n"
            << "  DENSITY, CC_DENSITY\\n"
            << "Phase 2 (Prefork+Fork+DivP1Late):\\n"
            << "  DIV_P1_PREFORK\\n"
            << "  Fork(AsyncWorker):\\n"
            << "    A: SET_BARO, VISC_BC, VEL_FLUX\\n"
            << "       PARTICLE_MOMENTUM\\n"
            << "    B: WALL_BC, DIV_P1_EARLY_B\\n"
            << "  DIV_P1_LATE_B\\n"
            << "Phase 3 (DivPart2):\\n"
            << "  DIV_P2_BLOCK_KERNEL";
        return oss.str();
    }
};

#endif // PRED_PREFORK_DIV_TASK_H
