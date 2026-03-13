#ifndef WALLBC_STATE_H
#define WALLBC_STATE_H

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
            double dt_bc = 0.0;
            int call_ht_1d = 0;

            // CALL_HT_1D only fires during the corrector phase (wall.f90:94)
            // WALL_COUNTER is only incremented before the corrector WALL_BC (main.f90:868)
            if (collected_[0]->phase == 1) { // corrector
                dt_bc = fds_compute_wall_bc_dt_bc(collected_[0]->t);
                fds_increment_wall_counter();
                call_ht_1d = fds_check_call_ht_1d();

                if (call_ht_1d) {
                    fds_update_bc_clock(collected_[0]->t);
                }
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
/// Gathers all N kernel results, runs sequential finalization, and emits
/// MeshData tokens downstream.
/// Finalization includes HAS_BACK_MESH processing, thin wall heat transfer, and particle off-gassing.
///
/// Meshes are placed directly at their correct position using NM as the index.
///
/// Flow: Collects N WallBCWork -> Sequential finalization -> Emits N MeshData
class WallBCCollector
    : public hh::AbstractState<1, WallBCWork, MeshData> {
public:
    explicit WallBCCollector(int nmeshes)
        : hh::AbstractState<1, WallBCWork, MeshData>(),
          nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<WallBCWork> work) override {
        collected_[work->nm - nmOffset_] = work;
        ++count_;

        if (count_ == nmeshes_) {
            // Sequential finalization: HAS_BACK_MESH cells, thin walls, particle off-gassing
            for (auto &w : collected_) {
                fds_wall_bc_finalize(w->nm, w->t, w->dt_bc, w->call_ht_1d);
            }
            // Reset WALL_COUNTER after WALL_BC loop (main.f90:872) — corrector only
            if (collected_[0]->originalMeshData->phase == 1) {
                fds_reset_wall_counter();
            }

            for (auto &w : collected_) {
                this->addResult(w->originalMeshData);
            }

            std::fill(collected_.begin(), collected_.end(), nullptr);
            count_ = 0;
        }
    }

private:
    int nmeshes_;
    int nmOffset_;
    int count_ = 0;
    std::vector<std::shared_ptr<WallBCWork>> collected_;
};

#endif // WALLBC_STATE_H
