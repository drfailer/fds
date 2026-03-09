#ifndef CORR_STEP1_STATE_H
#define CORR_STEP1_STATE_H

#include <algorithm>
#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/corr_step1_data.h"

/// Orchestrator state for corrector step 1 sub-graph.
///
/// Collects all N mesh tokens and dispatches parallel work tokens.
///
/// Flow: Collects N MeshData -> Emits N CorrStep1Work
class CorrStep1Orchestrator
    : public hh::AbstractState<1, MeshData, CorrStep1Work> {
public:
    explicit CorrStep1Orchestrator(int nmeshes)
        : hh::AbstractState<1, MeshData, CorrStep1Work>(),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
            for (auto &md : collected_) {
                auto work = std::make_shared<CorrStep1Work>(
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

/// Collector state for corrector step 1 sub-graph.
///
/// Gathers all N kernel results, sorts by mesh index for deterministic
/// ordering, and emits MeshData tokens downstream.
///
/// Flow: Collects N CorrStep1Work -> Emits N MeshData
class CorrStep1Collector
    : public hh::AbstractState<1, CorrStep1Work, MeshData> {
public:
    explicit CorrStep1Collector(int nmeshes)
        : hh::AbstractState<1, CorrStep1Work, MeshData>(),
          nmeshes_(nmeshes) {
        results_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<CorrStep1Work> work) override {
        results_.push_back(work);

        if (static_cast<int>(results_.size()) == nmeshes_) {
            // Sort by mesh index for deterministic ordering
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
    std::vector<std::shared_ptr<CorrStep1Work>> results_;
};

#endif // CORR_STEP1_STATE_H
