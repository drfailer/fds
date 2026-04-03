#ifndef VELOCITY_CORRECTOR_STATE_H
#define VELOCITY_CORRECTOR_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Orchestrator task for CC_IBM pre-processing before velocity corrector kernel.
///
/// Collects all N mesh tokens, then for each mesh runs:
///   1. CC_PROJECT_VELOCITY(STORE=TRUE) — store projected velocities
///   2. WALL_VELOCITY_NO_GRADH(STORE=TRUE) — store wall velocities for sparse solvers
///
/// Runs on a single thread.
class VelocityCorrectorCCOrchestrator
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit VelocityCorrectorCCOrchestrator(int nmeshes)
        : hh::AbstractTask<1, MeshData, MeshData>("VelCorrCCOrch", 1),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
            for (auto &md : collected_) {
                fds_cc_project_velocity_kernel(md->nm, md->dt, 1, 0);  // store=1, predictor=0
                fds_wall_velocity_no_gradh_kernel(md->nm, md->dt, 1, 0);  // store=1, predictor=0
            }

            for (auto &md : collected_) {
                this->addResult(md);
            }

            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

/// Collector task for CC_IBM post-processing after velocity corrector kernel.
///
/// Gathers all N kernel results, then for each mesh runs:
///   1. CC_PROJECT_VELOCITY(STORE=FALSE) — apply projected velocities
///   2. WALL_VELOCITY_NO_GRADH(STORE=FALSE) — restore wall velocities for sparse solvers
///
/// Runs on a single thread.
class VelocityCorrectorCCCollector
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit VelocityCorrectorCCCollector(int nmeshes)
        : hh::AbstractTask<1, MeshData, MeshData>("VelCorrCCCollector", 1),
          nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_[data->nm - nmOffset_] = data;
        ++count_;

        if (count_ == nmeshes_) {
            for (auto &md : collected_) {
                fds_cc_project_velocity_kernel(md->nm, md->dt, 0, 0);  // store=0, predictor=0
                fds_wall_velocity_no_gradh_kernel(md->nm, md->dt, 0, 0);  // store=0, predictor=0
            }

            for (auto &md : collected_) {
                this->addResult(md);
            }

            std::fill(collected_.begin(), collected_.end(), nullptr);
            count_ = 0;
        }
    }

private:
    int nmeshes_;
    int nmOffset_;
    int count_ = 0;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

#endif // VELOCITY_CORRECTOR_STATE_H
