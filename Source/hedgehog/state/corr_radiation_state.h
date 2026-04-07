#ifndef CORR_RADIATION_STATE_H
#define CORR_RADIATION_STATE_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../data/corr_radiation_data.h"
#include "../fds_fortran_interface.h"

#include <vector>

/// Orchestrator task: collects N MeshData tokens, dispatches parallel radiation work.
/// No sequential pre-processing needed (MESH_EXCHANGE(6) already done).
///
/// Runs on a single thread.
class CorrRadiationOrchestrator
    : public hh::AbstractTask<1, MeshData, CorrRadiationWork> {
public:
    explicit CorrRadiationOrchestrator(int nmeshes)
        : hh::AbstractTask<1, MeshData, CorrRadiationWork>("CorrRadOrch", 1),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nmeshes_) {
            for (auto &md : collected_) {
                this->bufferResult(std::make_shared<CorrRadiationWork>(
                    md->nm, md->t, 1, md));
            }
            this->flushResults<CorrRadiationWork>();
            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

/// Collector task: gathers N CorrRadiationWork results, accumulates global
/// RAD_Q_SUM/KFST4_SUM, emits N MeshData tokens downstream.
///
/// Uses direct indexed placement (nm - offset) to avoid sorting.
/// Runs on a single thread.
class CorrRadiationCollector
    : public hh::AbstractTask<1, CorrRadiationWork, MeshData> {
public:
    explicit CorrRadiationCollector(int nmeshes)
        : hh::AbstractTask<1, CorrRadiationWork, MeshData>("CorrRadCollector", 1),
          nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<CorrRadiationWork> work) override {
        collected_[work->nm - nmOffset_] = work;
        ++count_;
        if (count_ == nmeshes_) {
            // Accumulate per-mesh partial sums into global variables
            for (auto &w : collected_) {
                fds_accumulate_rad_sums(
                    w->radQSumPartial, w->kfst4SumPartial);
            }
            for (auto &w : collected_) {
                this->bufferResult(w->originalMeshData);
            }
            this->flushResults<MeshData>();
            std::fill(collected_.begin(), collected_.end(), nullptr);
            count_ = 0;
        }
    }

private:
    int nmeshes_;
    int nmOffset_;
    int count_ = 0;
    std::vector<std::shared_ptr<CorrRadiationWork>> collected_;
};

#endif // CORR_RADIATION_STATE_H
