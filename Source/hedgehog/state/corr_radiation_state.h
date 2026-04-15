#ifndef CORR_RADIATION_STATE_H
#define CORR_RADIATION_STATE_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/corr_radiation_data.h"
#include "../fds_fortran_interface.h"

#include <vector>

/// Collector task: gathers N CorrRadiationWork results, accumulates global
/// RAD_Q_SUM/KFST4_SUM, emits 1 BarrierData downstream.
///
/// Uses direct indexed placement (nm - offset) to avoid sorting.
/// Runs on a single thread.
class CorrRadiationCollector
    : public hh::AbstractTask<1, CorrRadiationWork, BarrierData> {
public:
    explicit CorrRadiationCollector(int nmeshes)
        : hh::AbstractTask<1, CorrRadiationWork, BarrierData>("CorrRadCollector", 1),
          nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<CorrRadiationWork> work) override {
        collected_[work->nm - nmOffset_] = work;
        if (++count_ == nmeshes_) {
            // Accumulate per-mesh partial sums into global variables
            for (auto &w : collected_) {
                fds_accumulate_rad_sums(
                    w->radQSumPartial, w->kfst4SumPartial);
            }
            auto bd = std::make_shared<BarrierData>();
            bd->meshes.resize(nmeshes_);
            for (auto &w : collected_) {
                bd->meshes[w->nm - nmOffset_] = w->originalMeshData;
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
