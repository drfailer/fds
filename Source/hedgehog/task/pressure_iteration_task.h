#ifndef PRESSURE_ITERATION_TASK_H
#define PRESSURE_ITERATION_TASK_H

#include <hedgehog/hedgehog.h>
#include <thread_utils/thread_pool.hpp>
#include <algorithm>
#include <sstream>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/termination_data.h"
#include "../fds_fortran_interface.h"

namespace pressure_iter_detail {

struct PressureBatchCtx {
    std::shared_ptr<MeshData<>> *meshes;
    int batchSize;
    int total;
    int presFlag;
    bool ccIBM;
    bool needsBaroclinic;
    int pressureIterations;
};

inline void baroclinicBatch(void *raw, TU_i64 batchIdx) {
    auto *ctx = static_cast<PressureBatchCtx *>(raw);
    int start = static_cast<int>(batchIdx) * ctx->batchSize;
    int end = std::min(start + ctx->batchSize, ctx->total);
    for (int i = start; i < end; ++i) {
        auto &md = ctx->meshes[i];
        if (ctx->needsBaroclinic) {
            fds_baroclinic_correction(md->t, md->nm);
        }
        if (ctx->ccIBM) {
            fds_cc_no_flux(md->dt, md->nm, 1);
            fds_cc_exchange_prepare_fn(md->nm);
        }
    }
}

inline void exchangeCopyBatch(void *raw, TU_i64 batchIdx) {
    auto *ctx = static_cast<PressureBatchCtx *>(raw);
    int start = static_cast<int>(batchIdx) * ctx->batchSize;
    int end = std::min(start + ctx->batchSize, ctx->total);
    for (int i = start; i < end; ++i) {
        int nm = ctx->meshes[i]->nm;
        int ncount = fds_flux_get_neighbor_count(nm);
        for (int j = 0; j < ncount; ++j) {
            int nom = fds_flux_get_neighbor_mesh(nm, j + 1);
            if (fds_flux_has_send_cells(nm, nom)) {
                fds_flux_copy_neighbor_ts(nm, nom);
            }
        }
    }
}

inline void solveBatch(void *raw, TU_i64 batchIdx) {
    auto *ctx = static_cast<PressureBatchCtx *>(raw);
    int start = static_cast<int>(batchIdx) * ctx->batchSize;
    int end = std::min(start + ctx->batchSize, ctx->total);
    for (int i = start; i < end; ++i) {
        auto &md = ctx->meshes[i];
        if (ctx->needsBaroclinic || ctx->pressureIterations == 1) {
            if (ctx->ccIBM) {
                fds_cc_match_velocity_flux(md->nm);
            } else {
                fds_match_velocity_flux_kernel(md->nm);
            }
        }
        fds_no_flux_kernel(md->nm, md->dt);
        if (ctx->ccIBM) {
            fds_cc_no_flux(md->dt, md->nm, 0);
        }
        if (ctx->pressureIterations == 1) {
            fds_pressure_iteration_zero_wall_work1(md->nm);
        }
        fds_pressure_solver_compute_rhs_kernel(md->nm, md->t, md->dt);
        if (ctx->presFlag == 3) {
            fds_ulmat_solver_kernel(md->nm, md->t, md->dt);
            fds_ulmat_check_residuals_kernel(md->nm);
        } else {
            fds_pressure_solver_fft_kernel(md->nm);
            fds_pressure_check_residuals_kernel(md->nm);
        }
    }
}

inline void velErrorBatch(void *raw, TU_i64 batchIdx) {
    auto *ctx = static_cast<PressureBatchCtx *>(raw);
    int start = static_cast<int>(batchIdx) * ctx->batchSize;
    int end = std::min(start + ctx->batchSize, ctx->total);
    for (int i = start; i < end; ++i) {
        auto &md = ctx->meshes[i];
        fds_compute_velocity_error_kernel(md->nm, md->dt);
        if (ctx->ccIBM) {
            fds_cc_compute_velocity_error(md->dt, md->nm);
        }
    }
}

} // namespace pressure_iter_detail

/// Single sequential task replacing the PressureIterationSubgraph.
///
/// Collects N MeshData tokens (from predictor or corrector), runs the
/// full pressure iteration loop internally, and emits converged results.
/// Parallel kernels (baroclinic, exchange copy, solve, velocity error)
/// use TU_ThreadPool with batched dispatch.
///
/// Shared between predictor and corrector via MeshData::phase routing.
/// TerminationData breaks the graph cycle for clean shutdown.
class PressureIterationTask : public hh::AbstractTask<3,
    MeshData<MeshState::PredictorPressure>,
    MeshData<MeshState::CorrectorPressure>,
    TerminationData,
    MeshData<MeshState::PredictorPressure>,
    MeshData<MeshState::CorrectorPressure>>
{
    using TaskBase = hh::AbstractTask<3,
        MeshData<MeshState::PredictorPressure>,
        MeshData<MeshState::CorrectorPressure>,
        TerminationData,
        MeshData<MeshState::PredictorPressure>,
        MeshData<MeshState::CorrectorPressure>>;

public:
    PressureIterationTask(int nmeshes, size_t threadCount, int presFlag)
        : TaskBase("PressureIteration", 1),
          nmeshes_(nmeshes),
          nmOffset_(fds_get_lower_mesh_index()),
          threadCount_(static_cast<int>(threadCount)),
          presFlag_(presFlag),
          ccIBM_(fds_is_cc_ibm() != 0) {
        collected_.resize(nmeshes, nullptr);
        int poolThreads = threadCount_ - 1;
        if (poolThreads > 0) {
            tu_tp_init(&pool_, static_cast<TU_u64>(poolThreads));
        }
    }

    ~PressureIterationTask() override {
        if (threadCount_ > 1) {
            tu_tp_fini(&pool_);
        }
    }

    void execute(std::shared_ptr<MeshData<MeshState::PredictorPressure>> md) override {
        collect(retag<MeshState::Default>(md));
    }

    void execute(std::shared_ptr<MeshData<MeshState::CorrectorPressure>> md) override {
        collect(retag<MeshState::Default>(md));
    }

    void execute(std::shared_ptr<TerminationData>) override {
        isDone_ = true;
    }

    [[nodiscard]] bool canTerminate() const override {
        return isDone_;
    }

private:
    void collect(std::shared_ptr<MeshData<>> md) {
        collected_[md->nm - nmOffset_] = md;
        if (++count_ < nmeshes_) return;
        count_ = 0;
        runPressureLoop();
    }

    void runPressureLoop() {
        for (;;) {
            int pIter = collected_[0]->pressure_iterations;
            bool iterBaro = collected_[0]->iterate_baroclinic;

            pressure_iter_detail::PressureBatchCtx ctx{
                collected_.data(),
                (nmeshes_ + threadCount_ - 1) / threadCount_,
                nmeshes_,
                presFlag_,
                ccIBM_,
                iterBaro,
                pIter};

            dispatchParallel(pressure_iter_detail::baroclinicBatch, ctx);
            dispatchParallel(pressure_iter_detail::exchangeCopyBatch, ctx);
            dispatchParallel(pressure_iter_detail::solveBatch, ctx);
            dispatchParallel(pressure_iter_detail::exchangeCopyBatch, ctx);
            dispatchParallel(pressure_iter_detail::velErrorBatch, ctx);

            int converged;
            if (fds_iterate_pressure()) {
                fds_set_pressure_iterations(pIter);
                fds_set_iterate_baroclinic_term(iterBaro ? 1 : 0);
                fds_pressure_iteration_check_convergence(
                    collected_[0]->t, collected_[0]->dt);
                converged = fds_pressure_iteration_converged();
                iterBaro = (fds_pressure_iteration_needs_baroclinic() != 0);
            } else {
                converged = 1;
            }

            if (converged) {
                emitConverged();
                return;
            }

            int nextPIter = pIter + 1;
            int totalPI = fds_get_total_pressure_iterations() + 1;
            fds_set_pressure_iterations(nextPIter);
            fds_set_iterate_baroclinic_term(iterBaro ? 1 : 0);
            fds_set_total_pressure_iterations(totalPI);

            for (auto &md : collected_) {
                md->pressure_iterations = nextPIter;
                md->iterate_baroclinic = iterBaro;
            }
        }
    }

    void emitConverged() {
        int phase = collected_[0]->phase;
        if (phase == 0) {
            fds_init_change_time_step(collected_[0]->dt);
            for (auto &md : collected_) {
                this->addResult(retag<MeshState::PredictorPressure>(md));
                md = nullptr;
            }
        } else {
            for (auto &md : collected_) {
                this->addResult(retag<MeshState::CorrectorPressure>(md));
                md = nullptr;
            }
        }
    }

    void dispatchParallel(tu_exec_func_t batchFn,
                          pressure_iter_detail::PressureBatchCtx &ctx) {
        int poolThreads = threadCount_ - 1;

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
            << "Pressure iteration loop:\\n"
            << "  Baroclinic (parallel)\\n"
            << "  MESH_EXCHANGE(5)\\n"
            << "  Solve (parallel)\\n"
            << "  VelError (parallel)\\n"
            << "  Convergence check";
        return oss.str();
    }

    int nmeshes_, nmOffset_, count_ = 0, threadCount_;
    int presFlag_;
    bool ccIBM_;
    bool isDone_ = false;
    TU_ThreadPool pool_{};
    std::vector<std::shared_ptr<MeshData<>>> collected_;
};

#endif // PRESSURE_ITERATION_TASK_H
