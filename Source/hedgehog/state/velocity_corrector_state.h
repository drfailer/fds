#ifndef VELOCITY_CORRECTOR_STATE_H
#define VELOCITY_CORRECTOR_STATE_H

#include <algorithm>
#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Orchestrator for CC_IBM pre-processing before velocity corrector kernel.
///
/// Collects all N mesh tokens, then for each mesh runs:
///   1. CC_PROJECT_VELOCITY(STORE=TRUE) — store projected velocities
///   2. WALL_VELOCITY_NO_GRADH(STORE=TRUE) — store wall velocities for sparse solvers
///
/// This matches velo.f90 VELOCITY_CORRECTOR (lines 646-656).
class VelocityCorrectorCCOrchestrator
    : public hh::AbstractState<1, MeshData, MeshData> {
public:
    explicit VelocityCorrectorCCOrchestrator(int nmeshes)
        : nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
            for (auto &md : collected_) {
                fds_cc_project_velocity(md->nm, md->dt, 1);  // STORE=.TRUE.
                fds_wall_velocity_no_gradh(md->nm, md->dt, 1);  // STORE=.TRUE.
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

/// Collector for CC_IBM post-processing after velocity corrector kernel.
///
/// Gathers all N kernel results, then for each mesh runs:
///   1. CC_PROJECT_VELOCITY(STORE=FALSE) — apply projected velocities
///   2. WALL_VELOCITY_NO_GRADH(STORE=FALSE) — restore wall velocities for sparse solvers
///
/// This matches velo.f90 VELOCITY_CORRECTOR (lines 660-670).
class VelocityCorrectorCCCollector
    : public hh::AbstractState<1, MeshData, MeshData> {
public:
    explicit VelocityCorrectorCCCollector(int nmeshes)
        : nmeshes_(nmeshes) {
        results_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        results_.push_back(data);

        if (static_cast<int>(results_.size()) == nmeshes_) {
            std::sort(results_.begin(), results_.end(),
                      [](const auto &a, const auto &b) { return a->nm < b->nm; });

            for (auto &md : results_) {
                fds_cc_project_velocity(md->nm, md->dt, 0);  // STORE=.FALSE.
                fds_wall_velocity_no_gradh(md->nm, md->dt, 0);  // STORE=.FALSE.
            }

            for (auto &md : results_) {
                this->addResult(md);
            }

            results_.clear();
            results_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> results_;
};

#endif // VELOCITY_CORRECTOR_STATE_H
