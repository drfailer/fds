#ifndef CORR_RADIATION_STATE_H
#define CORR_RADIATION_STATE_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/corr_radiation_data.h"
#include "../fds_fortran_interface.h"

#include <vector>

// Orchestrator: collects N MeshData tokens, dispatches parallel radiation work.
// No sequential pre-processing needed (MESH_EXCHANGE(6) already done).
class CorrRadiationOrchestrator
    : public hh::AbstractState<1, MeshData, CorrRadiationWork> {
public:
    explicit CorrRadiationOrchestrator(int nmeshes)
        : nmeshes_(nmeshes) { collected_.reserve(nmeshes); }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nmeshes_) {
            for (auto &md : collected_) {
                this->addResult(std::make_shared<CorrRadiationWork>(
                    md->nm, md->t, 1, md));
            }
            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

// Collector: gathers N CorrRadiationWork results, accumulates global
// RAD_Q_SUM/KFST4_SUM, sorts by nm, emits single BarrierData.
//
// Uses direct indexed placement (nm - offset) to avoid sorting.
class CorrRadiationCollector
    : public hh::AbstractState<1, CorrRadiationWork, BarrierData> {
public:
    explicit CorrRadiationCollector(int nmeshes)
        : nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()) {
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
            auto bd = std::make_shared<BarrierData>();
            bd->meshes.reserve(nmeshes_);
            for (auto &w : collected_) {
                bd->meshes.push_back(w->originalMeshData);
            }
            this->addResult(bd);
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
