#ifndef PRED_DIV_SEQUENTIAL_TASK_H
#define PRED_DIV_SEQUENTIAL_TASK_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Packed sequential task: ForkJoin + DivP1Late + DivExchange + DivP2Pre +
/// GlobalMatrix+PressureInit.
///
/// Single-threaded. Collects 2N tokens from fork branches (ForkJoin), then runs
/// the barrier-separated divergence pipeline sequentially:
///   ForkJoin(2N→N) → DivP1Late_loop → DivExchange → DivP2Pre_loop →
///   GlobalMatrix+PressureInit → emit N tokens for DivPart2 parallel kernel.
///
/// DivP1Late and DivP2Pre are LIGHT kernels (weight 1) — running them
/// sequentially in a loop avoids a Hedgehog cycle between sequential and
/// parallel tasks while still eliminating 5 graph nodes.
class PredDivSequentialTask : public hh::AbstractTask<1,
    MeshData<>,       // ForkJoin input (from fork branches)
    MeshData<>>       // → parallel: DivPart2 kernel
{
public:
    PredDivSequentialTask(int nmeshes, bool useParallelPressure)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>(
              "PredDivSequential", 1),
          nmeshes_(nmeshes),
          useParallelPressure_(useParallelPressure),
          nmOffset_(fds_get_lower_mesh_index()) {
        forkJoinCounts_.resize(nmeshes, 0);
        forkJoinReady_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData<>> data) override {
        int idx = data->nm - nmOffset_;
        forkJoinCounts_[idx]++;
        if (forkJoinCounts_[idx] == 2) {
            forkJoinCounts_[idx] = 0;
            forkJoinCompleted_++;
            forkJoinReady_.push_back(data);
            if (forkJoinCompleted_ == nmeshes_) {
                forkJoinCompleted_ = 0;
                runBarrierPipeline();
            }
        }
    }

private:
    /// Run the barrier-separated divergence pipeline on the collected meshes.
    void runBarrierPipeline() {
        // Phase 1: DivP1Late kernel (per-mesh, sequential loop)
        for (auto &md : forkJoinReady_) {
            fds_divergence_part_1_late_b(md->nm, md->t, md->dt);
        }

        // Phase 2: DivExchange barrier (global)
        fds_exchange_divergence_info();

        // Phase 3: DivP2 preprocessing kernel (per-mesh, sequential loop)
        for (auto &md : forkJoinReady_) {
            fds_divergence_part_2_preprocessing(md->nm, md->dt);
        }

        // Phase 4: GlobalMatrix + PressureInit barrier (global)
        fds_global_matrix_reassign(0);
        if (useParallelPressure_) {
            fds_pressure_iteration_init();
            fds_pressure_iteration_increment();
        }

        // Emit N tokens for DivPart2 parallel kernel
        for (auto &md : forkJoinReady_) {
            this->addResult(md);
        }
        forkJoinReady_.clear();
    }

    int nmeshes_, nmOffset_;
    bool useParallelPressure_;
    int forkJoinCompleted_ = 0;
    std::vector<int> forkJoinCounts_;
    std::vector<std::shared_ptr<MeshData<>>> forkJoinReady_;
};

#endif // PRED_DIV_SEQUENTIAL_TASK_H
