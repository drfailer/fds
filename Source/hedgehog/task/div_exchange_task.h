#ifndef DIV_EXCHANGE_TASK_H
#define DIV_EXCHANGE_TASK_H

#include <hedgehog/hedgehog.h>
#include <thread_utils/thread_pool.hpp>
#include <algorithm>
#include <mutex>
#include <sstream>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/termination_data.h"
#include "../fds_fortran_interface.h"

namespace div_exchange_detail {

struct BatchCtx {
    std::shared_ptr<MeshData<MeshState::DivExch>> *meshes;
    int batchSize;
    int total;
};

inline void divP2PreBatch(void *raw, TU_i64 batchIdx) {
    auto *ctx = static_cast<BatchCtx *>(raw);
    int start = static_cast<int>(batchIdx) * ctx->batchSize;
    int end = std::min(start + ctx->batchSize, ctx->total);
    for (int i = start; i < end; ++i) {
        fds_divergence_part_2_preprocessing(ctx->meshes[i]->nm, ctx->meshes[i]->dt);
    }
}

} // namespace div_exchange_detail

/// Merged DivExchange + DivP2Pre + GlobalMatrix task.
///
/// Replaces 2 barriers and the DivP2Pre phase from PredPreforkDivTask/CorrDivParallelTask.
/// Collects N MeshData<DivExch>, runs:
///   1. fds_exchange_divergence_info() (local reduction + MPI_ALLREDUCE)
///   2. fds_divergence_part_2_preprocessing() per mesh (parallel via TU_ThreadPool)
///   3. fds_global_matrix_reassign(0)
///   4. If parallel pressure: fds_pressure_iteration_init/increment
/// Then emits N MeshData<DivPart2>.
///
/// In a structural cycle with PredPreforkDivTask/CorrDivParallelTask — needs
/// TerminationData + canTerminate() for graph shutdown.
template<MeshState PressureTag = MeshState::Default>
class DivExchangeTask
    : public hh::AbstractTask<2,
        MeshData<MeshState::DivExch>, TerminationData,
        MeshData<MeshState::DivPart2>> {

    using TaskBase = hh::AbstractTask<2,
        MeshData<MeshState::DivExch>, TerminationData,
        MeshData<MeshState::DivPart2>>;

public:
    DivExchangeTask(int nmeshes, size_t threadCount)
        : TaskBase("DivExchange", 1),
          nmeshes_(nmeshes),
          nmOffset_(fds_get_lower_mesh_index()),
          threadCount_(static_cast<int>(threadCount)),
          useParallelPressure_(PressureTag != MeshState::Default) {
        collected_.resize(nmeshes, nullptr);
        int poolThreads = threadCount_ - 1;
        if (poolThreads > 0) {
            tu_tp_init(&pool_, static_cast<TU_u64>(poolThreads));
        }
    }

    ~DivExchangeTask() override {
        if (threadCount_ > 1) {
            tu_tp_fini(&pool_);
        }
    }

    DivExchangeTask(DivExchangeTask const &) = delete;
    DivExchangeTask &operator=(DivExchangeTask const &) = delete;

    void execute(std::shared_ptr<MeshData<MeshState::DivExch>> data) override {
        collected_[data->nm - nmOffset_] = data;
        if (++count_ < nmeshes_) return;
        count_ = 0;

        fds_exchange_divergence_info();

        dispatchParallel(div_exchange_detail::divP2PreBatch);

        fds_global_matrix_reassign(0);

        if (useParallelPressure_) {
            bool iterBaro = fds_get_baroclinic() != 0;
            int totalPI = fds_get_total_pressure_iterations() + 1;
            fds_set_pressure_iterations(1);
            fds_set_iterate_baroclinic_term(iterBaro ? 1 : 0);
            fds_set_total_pressure_iterations(totalPI);
            for (auto &md : collected_) {
                md->pressure_iterations = 1;
                md->iterate_baroclinic = iterBaro;
            }
        }

        for (auto &md : collected_) {
            this->addResult(retag<MeshState::DivPart2>(md));
            md = nullptr;
        }
    }

    void execute(std::shared_ptr<TerminationData>) override {
        std::lock_guard<std::mutex> lk(mtx_);
        done_ = true;
    }

    [[nodiscard]] bool canTerminate() const override {
        std::lock_guard<std::mutex> lk(mtx_);
        return done_;
    }

    [[nodiscard]] std::string extraPrintingInformation() const override {
        std::ostringstream oss;
        oss << "ThreadPool: " << threadCount_
            << " (" << (threadCount_ - 1) << " pool + 1 task)\\n"
            << "EXCH_DIV_INFO\\n"
            << "DivP2Pre (parallel)\\n"
            << "GLOBAL_MATRIX_REASSIGN";
        if (useParallelPressure_) {
            oss << "\\nPRES_INIT+INCR";
        }
        return oss.str();
    }

private:
    void dispatchParallel(tu_exec_func_t batchFn) {
        int poolThreads = threadCount_ - 1;
        int batchSize = (nmeshes_ + threadCount_ - 1) / threadCount_;

        div_exchange_detail::BatchCtx ctx{
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

    int nmeshes_, nmOffset_, count_ = 0, threadCount_;
    bool useParallelPressure_, done_ = false;
    mutable std::mutex mtx_;
    TU_ThreadPool pool_{};
    std::vector<std::shared_ptr<MeshData<MeshState::DivExch>>> collected_;
};

#endif // DIV_EXCHANGE_TASK_H
