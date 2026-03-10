#ifndef PRED_STEP1_STATE_H
#define PRED_STEP1_STATE_H

#include <hedgehog/hedgehog.h>
#include <algorithm>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/pred_step1_data.h"
#include "../fds_fortran_interface.h"

/// Orchestrator for PredStep1 sub-graph (Pattern B).
/// Collects N MeshData tokens, runs sequential INSERT_ALL_PARTICLES for each mesh,
/// then dispatches parallel work for COMPUTE_VISCOSITY + MASS_FINITE_DIFFERENCES kernels.
class PredStep1Orchestrator : public hh::AbstractState<1, MeshData, PredStep1Work> {
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> collected_;
public:
    explicit PredStep1Orchestrator(int nmeshes)
        : nmeshes_(nmeshes) { collected_.reserve(nmeshes); }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // Sequential pre-processing: INSERT_ALL_PARTICLES (cross-mesh, global state)
            for (auto &md : collected_) {
                fds_insert_particles(md->t, md->nm);
            }
            // Dispatch parallel kernel work
            for (auto &md : collected_) {
                this->addResult(std::make_shared<PredStep1Work>(md->nm, md->t, md->dt, md));
            }
            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }
};

/// Collector for PredStep1 sub-graph.
/// Collects N PredStep1Work results, sorts by mesh index, emits N MeshData.
class PredStep1Collector : public hh::AbstractState<1, PredStep1Work, MeshData> {
    int nmeshes_;
    std::vector<std::shared_ptr<PredStep1Work>> results_;
public:
    explicit PredStep1Collector(int nmeshes)
        : nmeshes_(nmeshes) { results_.reserve(nmeshes); }

    void execute(std::shared_ptr<PredStep1Work> work) override {
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

#endif // PRED_STEP1_STATE_H
