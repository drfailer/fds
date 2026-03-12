#ifndef CORR_RADIATION_STATE_H
#define CORR_RADIATION_STATE_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/corr_radiation_data.h"
#include "../fds_fortran_interface.h"

#include <algorithm>
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
class CorrRadiationCollector
    : public hh::AbstractState<1, CorrRadiationWork, BarrierData> {
public:
    explicit CorrRadiationCollector(int nmeshes)
        : nmeshes_(nmeshes) { results_.reserve(nmeshes); }

    void execute(std::shared_ptr<CorrRadiationWork> work) override {
        results_.push_back(work);
        if (static_cast<int>(results_.size()) == nmeshes_) {
            // Accumulate per-mesh partial sums into global variables
            for (auto &w : results_) {
                fds_accumulate_rad_sums(
                    w->radQSumPartial, w->kfst4SumPartial);
            }
            // Sort by mesh index for deterministic ordering
            std::sort(results_.begin(), results_.end(),
                      [](const auto &a, const auto &b) {
                          return a->nm < b->nm;
                      });
            auto bd = std::make_shared<BarrierData>();
            bd->meshes.reserve(nmeshes_);
            for (auto &w : results_) {
                bd->meshes.push_back(w->originalMeshData);
            }
            this->addResult(bd);
            results_.clear();
            results_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<CorrRadiationWork>> results_;
};

#endif // CORR_RADIATION_STATE_H
