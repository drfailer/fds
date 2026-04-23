#ifndef TIMESTEP_TASKS_H
#define TIMESTEP_TASKS_H

#include <hedgehog/hedgehog.h>
#include <thread_utils/thread_pool.hpp>
#include <thread_utils/async_worker.hpp>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <memory>
#include <sstream>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../fds_fortran_interface.h"

namespace timestep_detail {

struct GlobalDumpCtx {
    int *icyc;
    double t, dt;
};

inline void globalDumpFn(void *raw, TU_i64) {
    auto *ctx = static_cast<GlobalDumpCtx *>(raw);
    fds_set_diagnostics(*ctx->icyc, ctx->t, ctx->dt);
    fds_exchange_global_outputs(ctx->t, ctx->dt);
    fds_update_controls(ctx->t, ctx->dt);
    fds_dump_global_outputs(ctx->t, ctx->dt);
    fds_write_strings(ctx->t, ctx->dt);
    fds_write_diagnostics(ctx->t, ctx->dt);
}

struct MeshDumpCtx {
    std::shared_ptr<MeshData<>> *meshes;
    int batchSize;
    int total;
};

inline void meshDumpBatch(void *raw, TU_i64 batchIdx) {
    auto *ctx = static_cast<MeshDumpCtx *>(raw);
    int start = static_cast<int>(batchIdx) * ctx->batchSize;
    int end = std::min(start + ctx->batchSize, ctx->total);
    for (int i = start; i < end; ++i) {
        auto &md = ctx->meshes[i];
        fds_dump_mesh_outputs_ts(md->t, md->dt, md->nm);
    }
}

} // namespace timestep_detail

/// Merged timestep task: replaces DumpTask + TimestepState/Manager.
///
/// Collects N MeshData<> from corrector (or N MeshData<Init> on first iteration),
/// runs end-of-corrector work (RTE, reduce, dump) and timestep loop logic
/// (stop_check, adjust_dt), then cycles MeshData<> back to predictor.
///
/// Dump I/O uses TU_ThreadPool (per-mesh dump) + TU_AsyncWorker (global dump fork)
/// so global and per-mesh I/O overlap.
///
/// canTerminate() breaks the main cycle when the simulation ends.
class TimestepTask : public hh::AbstractTask<2,
    MeshData<>,
    MeshData<MeshState::Init>,
    MeshData<>,
    BarrierData>
{
    using TaskBase = hh::AbstractTask<2,
        MeshData<>,
        MeshData<MeshState::Init>,
        MeshData<>,
        BarrierData>;

public:
    TimestepTask(int nmeshes, size_t threadCount, double tEnd,
                 std::shared_ptr<int> icyc)
        : TaskBase("Timestep", 1),
          nmeshes_(nmeshes),
          nmOffset_(fds_get_lower_mesh_index()),
          threadCount_(static_cast<int>(threadCount)),
          tEnd_(tEnd),
          icyc_(std::move(icyc)) {
        collected_.resize(nmeshes, nullptr);
        int poolThreads = threadCount_ - 1;
        if (poolThreads > 0) {
            tu_tp_init(&pool_, static_cast<TU_u64>(poolThreads));
        }
        tu_aw_init(&asyncWorker_);
    }

    ~TimestepTask() override {
        if (threadCount_ > 1) {
            tu_tp_fini(&pool_);
        }
        tu_aw_fini(&asyncWorker_);
    }

    void execute(std::shared_ptr<MeshData<MeshState::Init>> data) override {
        initCollected_.push_back(data);
        if (static_cast<int>(initCollected_.size()) == nmeshes_) {
            for (auto &md : initCollected_) {
                this->addResult(retag<MeshState::Default>(md));
            }
            initCollected_.clear();
        }
    }

    void execute(std::shared_ptr<MeshData<>> data) override {
        collected_[data->nm - nmOffset_] = data;
        if (++count_ < nmeshes_) return;
        count_ = 0;
        runTimestepLogic();
    }

    [[nodiscard]] bool canTerminate() const override { return done_; }

    [[nodiscard]] std::string extraPrintingInformation() const override {
        std::ostringstream oss;
        oss << "ThreadPool: " << threadCount_
            << " (" << (threadCount_ - 1) << " pool + 1 task)\\n"
            << "RTE_SOURCE_CORRECTION\\n"
            << "REDUCE_HRR_MASS\\n"
            << "Dump: global(async) || mesh(pool)\\n"
            << "STOP_CHECK + cycle\\n"
            << std::fixed << std::setprecision(3) << totalTime_ << "s"
            << " / " << invocations_ << " calls";
        if (invocations_ > 0)
            oss << " / avg " << std::setprecision(3)
                << (totalTime_ * 1000.0 / invocations_) << "ms";
        if (dumpCount_ > 0)
            oss << "\\nPool dump: " << std::setprecision(3) << poolTime_ << "s"
                << " / " << dumpCount_ << " dumps"
                << " / avg " << std::setprecision(3)
                << (poolTime_ * 1000.0 / dumpCount_) << "ms";
        return oss.str();
    }

private:
    void runTimestepLogic() {
        auto tStart = std::chrono::steady_clock::now();

        double t = collected_[0]->t;
        double dt = collected_[0]->dt;

        fds_rte_source_correction();
        fds_reduce_hrr_mass(dt);

        bool anyDump = false;
        for (auto &md : collected_) {
            bool dump = false;
            fds_check_dump_schedule(md->t, md->nm, &dump);
            if (dump) { anyDump = true; break; }
        }

        if (anyDump) {
            globalCtx_ = {icyc_.get(), t, dt};
            tu_aw_exec(&asyncWorker_, timestep_detail::globalDumpFn,
                       &globalCtx_, 0);

            auto tPool = std::chrono::steady_clock::now();
            dispatchMeshDump();
            auto tPoolDone = std::chrono::steady_clock::now();
            poolTime_ += std::chrono::duration<double>(tPoolDone - tPool).count();

            tu_aw_wait(&asyncWorker_);
            ++dumpCount_;
        } else {
            timestep_detail::GlobalDumpCtx ctx{icyc_.get(), t, dt};
            timestep_detail::globalDumpFn(&ctx, 0);
        }

        fds_stop_check(1, t, dt);

        int stopStatus = fds_get_stop_status();
        if (t >= tEnd_ || stopStatus != 0) {
            done_ = true;
            auto bd = std::make_shared<BarrierData>();
            bd->done = true;
            this->addResult(bd);
            for (auto &md : collected_) { md = nullptr; }
            ++invocations_;
            totalTime_ += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - tStart).count();
            return;
        }

        fds_set_predictor(1);
        fds_set_first_pass(1);
        double newDt = fds_adjust_dt(t, dt);
        ++(*icyc_);
        fds_set_icyc(*icyc_);

        for (auto &md : collected_) {
            md->phase = 0;
            md->dt = newDt;
            md->firstPass = true;
            md->dt_bc = 0.0;
            md->call_ht_1d = 0;
        }

        this->batchAddResult(collected_);
        for (auto &md : collected_) { md = nullptr; }

        ++invocations_;
        totalTime_ += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - tStart).count();
    }

    void dispatchMeshDump() {
        int poolThreads = threadCount_ - 1;
        int batchSize = (nmeshes_ + threadCount_ - 1) / threadCount_;

        timestep_detail::MeshDumpCtx ctx{
            collected_.data(), batchSize, nmeshes_};

        if (poolThreads > 0) {
            std::vector<TU_ExecData> jobs(poolThreads);
            for (int i = 0; i < poolThreads; ++i) {
                jobs[i] = {timestep_detail::meshDumpBatch, &ctx,
                           static_cast<TU_i64>(i)};
            }
            TU_OperationHandle handle{};
            tu_tp_lauch(&pool_, jobs.data(),
                        static_cast<size_t>(poolThreads), &handle);

            timestep_detail::meshDumpBatch(&ctx,
                static_cast<TU_i64>(poolThreads));

            tu_tp_op_wait(&handle);
        } else {
            timestep_detail::meshDumpBatch(&ctx, 0);
        }
    }

    int nmeshes_, nmOffset_, count_ = 0, threadCount_;
    double tEnd_;
    std::shared_ptr<int> icyc_;
    bool done_ = false;
    TU_ThreadPool pool_{};
    TU_AsyncWorker asyncWorker_{};
    timestep_detail::GlobalDumpCtx globalCtx_{};
    std::vector<std::shared_ptr<MeshData<>>> collected_;
    std::vector<std::shared_ptr<MeshData<MeshState::Init>>> initCollected_;
    double totalTime_ = 0.0, poolTime_ = 0.0;
    int invocations_ = 0, dumpCount_ = 0;
};

#endif // TIMESTEP_TASKS_H
