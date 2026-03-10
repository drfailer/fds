#ifndef CORR_CONDENS_STATE_H
#define CORR_CONDENS_STATE_H

#include <hedgehog/hedgehog.h>
#include <algorithm>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/corr_condens_data.h"

/// Orchestrator for CorrCondens sub-graph (Pattern A — pure kernel).
/// Collects N MeshData tokens, dispatches parallel condensation kernel work.
class CorrCondensOrchestrator : public hh::AbstractState<1, MeshData, CorrCondensWork> {
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> collected_;
public:
    explicit CorrCondensOrchestrator(int nmeshes)
        : nmeshes_(nmeshes) { collected_.reserve(nmeshes); }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nmeshes_) {
            for (auto &md : collected_) {
                this->addResult(std::make_shared<CorrCondensWork>(md->nm, md->t, md->dt, md));
            }
            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }
};

/// Collector for CorrCondens sub-graph.
class CorrCondensCollector : public hh::AbstractState<1, CorrCondensWork, MeshData> {
    int nmeshes_;
    std::vector<std::shared_ptr<CorrCondensWork>> results_;
public:
    explicit CorrCondensCollector(int nmeshes)
        : nmeshes_(nmeshes) { results_.reserve(nmeshes); }

    void execute(std::shared_ptr<CorrCondensWork> work) override {
        results_.push_back(work);
        if (static_cast<int>(results_.size()) == nmeshes_) {
            std::sort(results_.begin(), results_.end(),
                      [](const auto &a, const auto &b) { return a->nm < b->nm; });
            for (auto &w : results_) {
                this->addResult(w->originalMeshData);
            }
            results_.clear();
            results_.reserve(nmeshes_);
        }
    }
};

#endif // CORR_CONDENS_STATE_H
