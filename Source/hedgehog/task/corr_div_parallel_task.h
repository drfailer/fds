#ifndef CORR_DIV_PARALLEL_TASK_H
#define CORR_DIV_PARALLEL_TASK_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Packed parallel task for the corrector divergence pipeline.
///
/// Two per-mesh parallel kernels in one multi-threaded task:
///   1. QRAddCopy:  MeshData<>        → QR addition + WORK1 copy   → MeshData<DivExch>
///   2. DivPart2:   MeshData<DivPart2>→ divergence part 2 block     → MeshData<PressureTag>
///
/// Between Phase 1 and Phase 2, DivExchangeTask collects N tokens, runs
/// exchange_divergence_info + parallel DivP2Pre + global_matrix_reassign,
/// then re-emits with DivPart2 tag.
template<MeshState PressureTag = MeshState::Default>
class CorrDivParallelTask : public hh::AbstractTask<2,
    MeshData<>,                        // from Join2+MeshExch2 → QRAddCopy
    MeshData<MeshState::DivPart2>,     // from DivExchangeTask → DivPart2 kernel
    MeshData<MeshState::DivExch>,      // → DivExchangeTask
    MeshData<PressureTag>>             // → downstream (pressure or VelCorr)
{
    using TaskBase = hh::AbstractTask<2,
        MeshData<>, MeshData<MeshState::DivPart2>,
        MeshData<MeshState::DivExch>,
        MeshData<PressureTag>>;

public:
    explicit CorrDivParallelTask(size_t numThreads)
        : TaskBase("CorrDivParallel", numThreads) {}

    /// Phase 1: QR addition + WORK1 copy (from Join2+MeshExch2)
    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_divergence_part_1_add_qr_b(data->nm);
        fds_copy_work1_b_to_work1(data->nm);
        this->addResult(retag<MeshState::DivExch>(data));
    }

    /// Phase 2: DivPart2 block kernel (from DivExchangeTask)
    void execute(std::shared_ptr<MeshData<MeshState::DivPart2>> data) override {
        int kbar = fds_get_kbar(data->nm);
        fds_divergence_part_2_block_kernel(data->nm, data->dt, 1, kbar);
        this->addResult(retag<PressureTag>(data));
    }

    std::shared_ptr<TaskBase> copy() override {
        return std::make_shared<CorrDivParallelTask<PressureTag>>(
            this->numberThreads());
    }
};

#endif // CORR_DIV_PARALLEL_TASK_H
