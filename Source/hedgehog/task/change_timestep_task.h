#ifndef CHANGE_TIMESTEP_TASK_H
#define CHANGE_TIMESTEP_TASK_H

#include <hedgehog/hedgehog.h>
#include <thread_utils/thread_pool.hpp>
#include <algorithm>
#include <sstream>
#include <vector>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

namespace change_timestep_detail {

struct BatchKernelCtx {
    std::shared_ptr<MeshData<>> *meshes;
    int batchSize;
    int total;
};

inline void momDivBatch(void *raw, TU_i64 batchIdx) {
    auto *ctx = static_cast<BatchKernelCtx *>(raw);
    int start = static_cast<int>(batchIdx) * ctx->batchSize;
    int end = std::min(start + ctx->batchSize, ctx->total);
    for (int i = start; i < end; ++i) {
        fds_particle_momentum_kernel(ctx->meshes[i]->nm, ctx->meshes[i]->dt);
        fds_divergence_part_1_kernel(ctx->meshes[i]->nm, ctx->meshes[i]->t,
                                     ctx->meshes[i]->dt);
    }
}

inline void divP2Batch(void *raw, TU_i64 batchIdx) {
    auto *ctx = static_cast<BatchKernelCtx *>(raw);
    int start = static_cast<int>(batchIdx) * ctx->batchSize;
    int end = std::min(start + ctx->batchSize, ctx->total);
    for (int i = start; i < end; ++i) {
        int kbar = fds_get_kbar(ctx->meshes[i]->nm);
        fds_divergence_part_2_block_kernel(ctx->meshes[i]->nm,
                                           ctx->meshes[i]->dt, 1, kbar);
    }
}

inline void velPredBatch(void *raw, TU_i64 batchIdx) {
    auto *ctx = static_cast<BatchKernelCtx *>(raw);
    int start = static_cast<int>(batchIdx) * ctx->batchSize;
    int end = std::min(start + ctx->batchSize, ctx->total);
    for (int i = start; i < end; ++i) {
        auto &md = ctx->meshes[i];
        fds_velocity_predictor_kernel_only(md->nm, md->dt);
        fds_cc_project_velocity_kernel(md->nm, md->dt, 0, 1);
        fds_wall_velocity_no_gradh_kernel(md->nm, md->dt, 0, 1);
        fds_check_stability_kernel_only(md->nm, md->t + md->dt, md->dt);
    }
}

} // namespace change_timestep_detail

/// Single sequential task replacing the ChangeTimeStepSubgraph.
///
/// Collects N MeshData<> from VelocityPredictor, checks CFL condition,
/// and retries with smaller dt if needed. Parallel kernels (MomDiv, DivP2,
/// VelPred) use TU_ThreadPool with batched dispatch: pool threads + task
/// thread = threadCount total parallelism.
class ChangeTimeStepTask
    : public hh::AbstractTask<1, MeshData<MeshState::PostVelPred>, MeshData<MeshState::MeshExch3>> {
public:
    ChangeTimeStepTask(int nmeshes, size_t threadCount, bool ccIBM)
        : hh::AbstractTask<1, MeshData<MeshState::PostVelPred>, MeshData<MeshState::MeshExch3>>("ChangeTimeStep", 1),
          nmeshes_(nmeshes),
          nmOffset_(fds_get_lower_mesh_index()),
          threadCount_(static_cast<int>(threadCount)),
          ccIBM_(ccIBM) {
        collected_.resize(nmeshes, nullptr);
        int poolThreads = threadCount_ - 1;
        if (poolThreads > 0) {
            tu_tp_init(&pool_, static_cast<TU_u64>(poolThreads));
        }
    }

    ~ChangeTimeStepTask() override {
        if (threadCount_ > 1) {
            tu_tp_fini(&pool_);
        }
    }

    void execute(std::shared_ptr<MeshData<MeshState::PostVelPred>> tagged) override {
        auto data = retag<MeshState::Default>(tagged);
        collected_[data->nm - nmOffset_] = data;
        if (++count_ < nmeshes_) return;
        count_ = 0;

        fds_stop_check_zero();

        int needRetry = 0;
        double newDt = 0.0;
        fds_check_change_time_step(&needRetry, &newDt);

        if (!needRetry) {
            finishAndEmit();
            return;
        }

        double t = collected_[0]->t;

        for (;;) {
            doPreKernelWork(t, newDt);

            dispatchParallel(change_timestep_detail::momDivBatch);

            fds_exchange_divergence_info();
            // TODO: move this to parallel dispatch
            for (auto &md : collected_) {
                fds_divergence_part_2_preprocessing(md->nm, md->dt);
            }

            dispatchParallel(change_timestep_detail::divP2Batch);

            fds_pressure_iteration(collected_[0]->t, collected_[0]->dt);
            fds_init_change_time_step(collected_[0]->dt);

            dispatchParallel(change_timestep_detail::velPredBatch);

            fds_stop_check_zero();
            if (fds_get_stop_status() != 0) {
                finishAndEmit();
                return;
            }

            fds_check_change_time_step(&needRetry, &newDt);
            if (!needRetry) {
                finishAndEmit();
                return;
            }
        }
    }

private:
    void doPreKernelWork(double t, double dt) {
        fds_set_first_pass(0);

        for (auto &md : collected_) {
            md->dt = dt;
            md->firstPass = false;
        }

        for (auto &md : collected_) {
            fds_cc_restore_uvw_unlinked(md->nm);
            fds_density(t, dt, md->nm);
            fds_cc_density_ts(md->nm, t, dt);
        }

        fds_mesh_exchange(1);

        for (auto &md : collected_) {
            fds_set_baroclinic_false(md->nm);
            fds_viscosity_bc(md->nm, 0);
            fds_velocity_flux(t, dt, md->nm, 0);
        }

        fds_hvac_calc(t, dt, 0);
        fds_initialize_divergence_integrals();

        for (auto &md : collected_) {
            fds_wall_bc(t, dt, md->nm);
        }
    }

    void finishAndEmit() {
        if (ccIBM_) {
            fds_cc_end_step(collected_[0]->t, collected_[0]->dt, 0);
            fds_mesh_cc_exchange(3);
        }
        for (auto &md : collected_) {
            this->addResult(retag<MeshState::MeshExch3>(md));
            md = nullptr;
        }
    }

    void dispatchParallel(tu_exec_func_t batchFn) {
        int poolThreads = threadCount_ - 1;
        int batchSize = (nmeshes_ + threadCount_ - 1) / threadCount_;

        change_timestep_detail::BatchKernelCtx ctx{
            collected_.data(), batchSize, nmeshes_};

        if (poolThreads > 0) {
            std::vector<TU_ExecData> jobs(poolThreads);
            for (int i = 0; i < poolThreads; ++i) {
                jobs[i] = {batchFn, &ctx, static_cast<TU_i64>(i)};
            }
            TU_OperationHandle handle{};
            tu_tp_lauch(&pool_, jobs.data(),
                        static_cast<size_t>(poolThreads), &handle);

            batchFn(&ctx, static_cast<TU_i64>(poolThreads));

            tu_tp_op_wait(&handle);
        } else {
            batchFn(&ctx, 0);
        }
    }

    [[nodiscard]] std::string extraPrintingInformation() const override {
        std::ostringstream oss;
        oss << "ThreadPool: " << threadCount_
            << " (" << (threadCount_ - 1) << " pool + 1 task)\\n"
            << "CFL retry loop:\\n"
            << "  MomDiv (parallel)\\n"
            << "  EXCH_DIV_INFO\\n"
            << "  DivP2Pre (sequential)\\n"
            << "  DivP2Block (parallel)\\n"
            << "  PRESSURE_ITERATION\\n"
            << "  VelPred (parallel)";
        return oss.str();
    }

    int nmeshes_, nmOffset_, count_ = 0, threadCount_;
    bool ccIBM_;
    TU_ThreadPool pool_{};
    std::vector<std::shared_ptr<MeshData<>>> collected_;
};

#endif // CHANGE_TIMESTEP_TASK_H
