#ifndef WALLBC_STATE_H
#define WALLBC_STATE_H

#include <algorithm>
#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/wallbc_data.h"
#include "../fds_fortran_interface.h"

/// Orchestrator state for WallBC sub-graph.
///
/// Collects all N mesh tokens, computes global DT_BC/CALL_HT_1D state,
/// then dispatches parallel work. Per-mesh preprocessing (ASSIGN_GHOST_VALUE_KERNEL,
/// NEAR_SURFACE_GAS_VARIABLES, HEAT_TRANSFER_COEFFICIENT) has been moved to the
/// parallel kernel task.
///
/// Flow: Collects N MeshData -> Compute global state -> Emits N WallBCWork
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
            // Compute dt_bc and call_ht_1d from global Fortran state (same for all meshes)
            double dt_bc = fds_compute_wall_bc_dt_bc(collected_[0]->t);
            int call_ht_1d = fds_check_call_ht_1d();

            // If calling 1-D heat transfer, update BC_CLOCK
            if (call_ht_1d) {
                fds_update_bc_clock(collected_[0]->t);
            }

            // Dispatch parallel work (preprocessing + cell processing in kernel task)
            for (auto &md : collected_) {
                auto work = std::make_shared<WallBCWork>(
                    md->nm, md->t, md->dt, dt_bc, call_ht_1d, md);
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

/// Collector state for WallBC sub-graph (Pattern B - Sequential Post-Processing).
///
/// Gathers all N kernel results, runs sequential finalization, sorts by mesh index,
/// and emits MeshData tokens downstream.
/// Finalization includes HAS_BACK_MESH processing, thin wall heat transfer, and particle off-gassing.
///
/// Flow: Collects N WallBCWork -> Sequential finalization -> Emits N MeshData
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

            // Sequential finalization: HAS_BACK_MESH cells, thin walls, particle off-gassing
            for (auto &w : results_) {
                fds_wall_bc_finalize(w->nm, w->t, w->dt_bc, w->call_ht_1d);
            }

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
