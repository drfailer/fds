#ifndef DENSITY_PRED_STATE_H
#define DENSITY_PRED_STATE_H

#include <algorithm>
#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/density_pred_data.h"

/// Orchestrator state for density predictor sub-graph.
class DensityPredOrchestrator
    : public hh::AbstractState<1, MeshData, DensityPredWork> {
public:
    explicit DensityPredOrchestrator(int nmeshes)
        : hh::AbstractState<1, MeshData, DensityPredWork>(),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
            for (auto &md : collected_) {
                auto work = std::make_shared<DensityPredWork>(
                    md->nm, md->t, md->dt, md);
                this->addResult(work);
            }

            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

/// Collector state for density predictor sub-graph.
class DensityPredCollector
    : public hh::AbstractState<1, DensityPredWork, MeshData> {
public:
    explicit DensityPredCollector(int nmeshes)
        : hh::AbstractState<1, DensityPredWork, MeshData>(),
          nmeshes_(nmeshes) {
        results_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<DensityPredWork> work) override {
        results_.push_back(work);

        if (static_cast<int>(results_.size()) == nmeshes_) {
            std::sort(results_.begin(), results_.end(),
                      [](const auto &a, const auto &b) {
                          return a->nm < b->nm;
                      });

            for (auto &w : results_) {
                this->addResult(w->originalMeshData);
            }

            results_.clear();
            results_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<DensityPredWork>> results_;
};

#endif // DENSITY_PRED_STATE_H
