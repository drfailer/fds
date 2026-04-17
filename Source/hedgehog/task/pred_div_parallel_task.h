#ifndef PRED_DIV_PARALLEL_TASK_H
#define PRED_DIV_PARALLEL_TASK_H

#include <hedgehog/hedgehog.h>
#include <atomic>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/termination_data.h"
#include "../fds_fortran_interface.h"

/// Packed parallel "thread pool" task for the predictor divergence pipeline.
///
/// Merges 3 per-mesh parallel kernels into one multi-threaded task:
///   1. DivP1Late:  MeshData<>        → fds_divergence_part_1_late_b  → MeshData<DivExch>
///   2. DivP2Pre:   MeshData<DivP2Pre>→ fds_divergence_part_2_preprocessing → MeshData<GlobalMat>
///   3. DivPart2:   MeshData<DivPart2>→ fds_divergence_part_2_block_kernel  → MeshData<PressureTag>
///
/// Between these phases, sequential barriers (DivExchange, GlobalMatrix) collect
/// N tokens, run global operations, and re-emit with the next phase's tag.
/// Although all input/output types are distinct (no data-level cycle), Hedgehog
/// detects structural cycles at the node level. TerminationData + canTerminate()
/// break this at shutdown.
///
/// One thread pool serves all 3 phases, saving 2×solo(1) threads vs. separate tasks.
template<MeshState PressureTag = MeshState::Default>
class PredDivParallelTask : public hh::AbstractTask<4,
    MeshData<>,                        // from ForkJoin → DivP1Late kernel
    MeshData<MeshState::DivP2Pre>,     // from DivExchange → DivP2Pre kernel
    MeshData<MeshState::DivPart2>,     // from GlobalMatrix → DivPart2 kernel
    TerminationData,                   // termination signal
    MeshData<MeshState::DivExch>,      // → DivExchange barrier
    MeshData<MeshState::GlobalMat>,    // → GlobalMatrix barrier
    MeshData<PressureTag>>             // → downstream (pressure or VelPred)
{
    using TaskBase = hh::AbstractTask<4,
        MeshData<>, MeshData<MeshState::DivP2Pre>, MeshData<MeshState::DivPart2>,
        TerminationData,
        MeshData<MeshState::DivExch>, MeshData<MeshState::GlobalMat>,
        MeshData<PressureTag>>;

public:
    explicit PredDivParallelTask(size_t numThreads)
        : TaskBase("PredDivParallel", numThreads),
          done_(std::make_shared<std::atomic<bool>>(false)) {}

    /// Phase 1: DivP1Late kernel (from ForkJoin)
    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_divergence_part_1_late_b(data->nm, data->t, data->dt);
        this->addResult(retag<MeshState::DivExch>(data));
    }

    /// Phase 2: DivP2 preprocessing kernel (from DivExchange barrier)
    void execute(std::shared_ptr<MeshData<MeshState::DivP2Pre>> data) override {
        fds_divergence_part_2_preprocessing(data->nm, data->dt);
        this->addResult(retag<MeshState::GlobalMat>(data));
    }

    /// Phase 3: DivPart2 block kernel (from GlobalMatrix barrier)
    void execute(std::shared_ptr<MeshData<MeshState::DivPart2>> data) override {
        int kbar = fds_get_kbar(data->nm);
        fds_divergence_part_2_block_kernel(data->nm, data->dt, 1, kbar);
        this->addResult(retag<PressureTag>(data));
    }

    /// Termination signal — sets done flag so canTerminate() returns true.
    void execute(std::shared_ptr<TerminationData>) override {
        done_->store(true);
    }

    /// Break structural cycle at shutdown. During normal operation, the outer
    /// timestep loop keeps all predecessors alive. At shutdown, TerminationData
    /// sets done_=true, this node terminates, and the barrier nodes cascade.
    [[nodiscard]] bool canTerminate() const override {
        return done_->load();
    }

    std::shared_ptr<TaskBase> copy() override {
        auto c = std::make_shared<PredDivParallelTask<PressureTag>>(
            this->numberThreads());
        c->done_ = done_;  // share atomic flag across clones
        return c;
    }

    std::shared_ptr<std::atomic<bool>> done_;
};

#endif // PRED_DIV_PARALLEL_TASK_H
