#ifndef DIVERGENCE_PART2_STATE_H
#define DIVERGENCE_PART2_STATE_H

#include <algorithm>
#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/divergence_part2_data.h"

/// Orchestrator state for divergence part 2 sub-graph.
///
/// Collects all N mesh tokens and dispatches parallel work tokens.
///
/// Flow: Collects N MeshData -> Emits N DivergencePart2Work
class DivergencePart2Orchestrator
    : public hh::AbstractState<1, MeshData, DivergencePart2Work> {
public:
    explicit DivergencePart2Orchestrator(int nmeshes)
        : hh::AbstractState<1, MeshData, DivergencePart2Work>(),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
            for (auto &md : collected_) {
                auto work = std::make_shared<DivergencePart2Work>(
                    md->nm, md->dt, md);
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

/// Collector state for divergence part 2 sub-graph.
///
/// Gathers all N kernel results, sorts by mesh index for deterministic
/// ordering, and emits MeshData tokens downstream.
///
/// Flow: Collects N DivergencePart2Work -> Emits N MeshData
class DivergencePart2Collector
    : public hh::AbstractState<1, DivergencePart2Work, MeshData> {
public:
    explicit DivergencePart2Collector(int nmeshes)
        : hh::AbstractState<1, DivergencePart2Work, MeshData>(),
          nmeshes_(nmeshes) {
        results_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<DivergencePart2Work> work) override {
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
    std::vector<std::shared_ptr<DivergencePart2Work>> results_;
};

#endif // DIVERGENCE_PART2_STATE_H
