#ifndef VELOCITY_CORRECTOR_STATE_H
#define VELOCITY_CORRECTOR_STATE_H

#include <algorithm>
#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/velocity_corrector_data.h"

/// Orchestrator state for velocity corrector sub-graph.
///
/// This state collects all N mesh tokens at the start of the velocity
/// corrector sub-graph and dispatches parallel work tokens for kernel
/// execution. Sequential pre-processing (e.g., CC_IBM setup) can be
/// performed here before dispatching.
///
/// Flow: Collects N MeshData → Emits N VelocityCorrectorWork
class VelocityCorrectorOrchestrator
    : public hh::AbstractState<1, MeshData, VelocityCorrectorWork> {
public:
    explicit VelocityCorrectorOrchestrator(int nmeshes)
        : hh::AbstractState<1, MeshData, VelocityCorrectorWork>(),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // All mesh tokens collected - ready to dispatch parallel work

            // Sequential pre-processing would go here
            // (Currently none needed for velocity corrector)

            // Emit work tokens for parallel kernel execution
            for (auto &md : collected_) {
                auto work = std::make_shared<VelocityCorrectorWork>(
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

/// Collector state for velocity corrector sub-graph.
///
/// This state gathers all N kernel execution results and performs any
/// sequential post-processing (e.g., global diagnostics, reductions)
/// before emitting MeshData tokens to continue the graph flow.
///
/// Flow: Collects N VelocityCorrectorWork → Emits N MeshData
class VelocityCorrectorCollector
    : public hh::AbstractState<1, VelocityCorrectorWork, MeshData> {
public:
    explicit VelocityCorrectorCollector(int nmeshes)
        : hh::AbstractState<1, VelocityCorrectorWork, MeshData>(),
          nmeshes_(nmeshes) {
        results_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<VelocityCorrectorWork> work) override {
        results_.push_back(work);

        if (static_cast<int>(results_.size()) == nmeshes_) {
            // Sort by mesh index to guarantee deterministic ordering
            std::sort(results_.begin(), results_.end(),
                      [](const auto &a, const auto &b) { return a->nm < b->nm; });

            // Emit original MeshData tokens to continue graph flow
            for (auto &w : results_) {
                this->addResult(w->originalMeshData);
            }

            results_.clear();
            results_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<VelocityCorrectorWork>> results_;
};

#endif // VELOCITY_CORRECTOR_STATE_H
