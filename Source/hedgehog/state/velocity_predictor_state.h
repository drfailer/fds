#ifndef VELOCITY_PREDICTOR_STATE_H
#define VELOCITY_PREDICTOR_STATE_H

#include <algorithm>
#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/velocity_predictor_data.h"

/// Orchestrator state for velocity predictor sub-graph.
///
/// Collects all N mesh tokens and dispatches parallel work tokens for kernel
/// execution. Sequential pre-processing (e.g., CC_IBM setup) can be
/// performed here before dispatching.
///
/// Flow: Collects N MeshData -> Emits N VelocityPredictorWork
class VelocityPredictorOrchestrator
    : public hh::AbstractState<1, MeshData, VelocityPredictorWork> {
public:
    explicit VelocityPredictorOrchestrator(int nmeshes)
        : hh::AbstractState<1, MeshData, VelocityPredictorWork>(),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // All mesh tokens collected - ready to dispatch parallel work

            // Sequential pre-processing would go here
            // (e.g., CC_IBM CC_PROJECT_VELOCITY, WALL_VELOCITY_NO_GRADH
            //  if needed in the future)

            // Emit work tokens for parallel kernel execution
            for (auto &md : collected_) {
                auto work = std::make_shared<VelocityPredictorWork>(
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

/// Collector state for velocity predictor sub-graph.
///
/// Gathers all N kernel execution results and performs any sequential
/// post-processing before emitting MeshData tokens to continue the graph flow.
///
/// Flow: Collects N VelocityPredictorWork -> Emits N MeshData
class VelocityPredictorCollector
    : public hh::AbstractState<1, VelocityPredictorWork, MeshData> {
public:
    explicit VelocityPredictorCollector(int nmeshes)
        : hh::AbstractState<1, VelocityPredictorWork, MeshData>(),
          nmeshes_(nmeshes) {
        results_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<VelocityPredictorWork> work) override {
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
    std::vector<std::shared_ptr<VelocityPredictorWork>> results_;
};

#endif // VELOCITY_PREDICTOR_STATE_H
