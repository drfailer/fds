#ifndef PRED_PREFORK_DIV_TASK_H
#define PRED_PREFORK_DIV_TASK_H

#include <hedgehog/hedgehog.h>
#include <sstream>
#include <thread_utils/async_worker.hpp>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Merged predictor task: DivP1Prefork + Fork(DivSetup+PartMom || WallBC+DivP1Early).
///
/// Replaces 3 separate tasks (DivP1PreforkKernel, PredDivSetupPartMom,
/// PredWallBCDivEarly) + ForkJoinTask by using an AsyncWorker for the fork.
///
/// Per-mesh execution:
///   1. DIV_P1_PREFORK (main thread, sequential)
///   2. Fork via AsyncWorker:
///      - Worker thread: WALL_BC + DIV_P1_EARLY_B
///      - Main thread:   DivSetup + PARTICLE_MOMENTUM
///   3. Join: tu_aw_wait
///
/// Each copy() creates its own AsyncWorker — no sharing between Hedgehog
/// threads.  Real OS thread count = 2 x numThreads (main + worker per thread).
class PredPreforkDivTask
    : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
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
        : hh::AbstractTask<1, MeshData<>, MeshData<>>(
              "PredPreforkDiv", numThreads) {
        tu_aw_init(&worker_);
    }

    ~PredPreforkDivTask() override { tu_aw_fini(&worker_); }

    PredPreforkDivTask(PredPreforkDivTask const &) = delete;
    PredPreforkDivTask &operator=(PredPreforkDivTask const &) = delete;

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_divergence_part_1_prefork(data->nm, data->t, data->dt);

        tu_aw_exec(&worker_, wallBCDivEarlyWork, data.get(), data->nm);

        fds_set_baroclinic_false(data->nm);
        fds_viscosity_bc_kernel(data->nm, data->phase);
        fds_cc_velocity_bc_ts(data->t, data->nm, data->phase, 0);
        fds_velocity_flux_kernel(data->nm, data->t, data->dt, data->phase);
        fds_particle_momentum_kernel(data->nm, data->dt);

        tu_aw_wait(&worker_);

        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>>
    copy() override {
        return std::make_shared<PredPreforkDivTask>(this->numberThreads());
    }

    [[nodiscard]] std::string extraPrintingInformation() const override {
        std::ostringstream oss;
        oss << "DIV_P1_PREFORK\\n"
            << "Fork(AsyncWorker):\\n"
            << "  A: SET_BARO_FALSE, VISC_BC\\n"
            << "     CC_VEL_BC_TS, VEL_FLUX\\n"
            << "     PARTICLE_MOMENTUM\\n"
            << "  B: WALL_BC_PREPROC\\n"
            << "     WALL_BC_CELLS\\n"
            << "     WALL_BC_FINALIZE\\n"
            << "     DIV_P1_EARLY_B";
        return oss.str();
    }
};

#endif // PRED_PREFORK_DIV_TASK_H
