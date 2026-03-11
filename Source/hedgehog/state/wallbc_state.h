#ifndef WALLBC_STATE_H
#define WALLBC_STATE_H

#include <algorithm>
#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/wallbc_data.h"

/// Orchestrator state for WallBC sub-graph (Pattern A - Pure Kernel).
///
/// Collects all N mesh tokens and dispatches parallel work tokens.
/// Note: Preprocessing (ASSIGN_GHOST_VALUE, NEAR_SURFACE_GAS_VARIABLES, etc.)
/// must be handled before this sub-graph is invoked.
///
/// Flow: Collects N MeshData -> Emits N WallBCWork
class WallBCOrchestrator
    : public hh::AbstractState<1, MeshData, WallBCWork> {
public:
    explicit WallBCOrchestrator(int nmeshes)
        : hh::AbstractState<1, MeshData, WallBCWork>(),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // Dispatch parallel work for WALL_BC_PROCESS_CELLS_KERNEL
            // DT_BC and CALL_HT_1D are computed from global state (same for all meshes)
            // For now, we pass them as part of the work token
            // TODO: Extract these from Fortran global state or pass via MeshData
            for (auto &md : collected_) {
                auto work = std::make_shared<WallBCWork>(
                    md->nm, md->t, md->dt,
                    md->dt,  // dt_bc (placeholder - should be computed properly)
                    1,        // call_ht_1d (placeholder - should be computed properly)
                    md);
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

/// Collector state for WallBC sub-graph.
///
/// Gathers all N kernel results, sorts by mesh index for deterministic
/// ordering, and emits MeshData tokens downstream.
/// Note: Finalization (WALL_BC_FINALIZE) must be handled after this sub-graph.
///
/// Flow: Collects N WallBCWork -> Emits N MeshData
class WallBCCollector
    : public hh::AbstractState<1, WallBCWork, MeshData> {
public:
    explicit WallBCCollector(int nmeshes)
        : hh::AbstractState<1, WallBCWork, MeshData>(),
          nmeshes_(nmeshes) {
        results_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<WallBCWork> work) override {
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
    std::vector<std::shared_ptr<WallBCWork>> results_;
};

#endif // WALLBC_STATE_H
