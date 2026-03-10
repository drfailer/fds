#ifndef PRED_WALL_DIV_STATE_H
#define PRED_WALL_DIV_STATE_H

#include <hedgehog/hedgehog.h>
#include <algorithm>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/pred_wall_div_data.h"
#include "../fds_fortran_interface.h"

/// Orchestrator for PredWallDiv sub-graph (Pattern B).
/// Collects N MeshData tokens, runs sequential WALL_BC for each mesh (cross-mesh OMESH access),
/// then dispatches parallel work for PARTICLE_MOMENTUM + DIVERGENCE_PART_1 kernels.
class PredWallDivOrchestrator : public hh::AbstractState<1, MeshData, PredWallDivWork> {
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> collected_;
public:
    explicit PredWallDivOrchestrator(int nmeshes)
        : nmeshes_(nmeshes) { collected_.reserve(nmeshes); }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // Sequential pre-processing: WALL_BC (reads OMESH for ghost cells)
            for (auto &md : collected_) {
                fds_wall_bc(md->t, md->dt, md->nm);
            }
            // Dispatch parallel kernel work
            for (auto &md : collected_) {
                this->addResult(std::make_shared<PredWallDivWork>(md->nm, md->t, md->dt, md));
            }
            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }
};

/// Collector for PredWallDiv sub-graph.
class PredWallDivCollector : public hh::AbstractState<1, PredWallDivWork, MeshData> {
    int nmeshes_;
    std::vector<std::shared_ptr<PredWallDivWork>> results_;
public:
    explicit PredWallDivCollector(int nmeshes)
        : nmeshes_(nmeshes) { results_.reserve(nmeshes); }

    void execute(std::shared_ptr<PredWallDivWork> work) override {
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

#endif // PRED_WALL_DIV_STATE_H
