#ifndef PRED_WALL_DIV_STATE_H
#define PRED_WALL_DIV_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Orchestrator for PredWallDiv sub-graph (Pattern B).
/// Collects N MeshData tokens, runs sequential WALL_BC for each mesh (cross-mesh OMESH access),
/// then dispatches MeshData for parallel PARTICLE_MOMENTUM + DIVERGENCE_PART_1 kernels.
class PredWallDivOrchestrator : public hh::AbstractState<1, MeshData, MeshData> {
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
            // Dispatch for parallel kernel execution
            for (auto &md : collected_) {
                this->addResult(md);
            }
            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }
};

#endif // PRED_WALL_DIV_STATE_H
