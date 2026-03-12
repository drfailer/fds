#ifndef VELOCITY_BC_STATE_H
#define VELOCITY_BC_STATE_H

#include <algorithm>
#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../fds_fortran_interface.h"

/// Orchestrator for PredFinal sub-graph.
///
/// Collects N mesh tokens, runs sequential SYNTHETIC_TURBULENCE_IF_ENABLED
/// (SEM inflow — uses RANDOM_NUMBER, kept sequential), then dispatches
/// N MeshData tokens for parallel kernel processing.
class PredFinalOrchestrator
    : public hh::AbstractState<1, MeshData, MeshData> {
public:
    explicit PredFinalOrchestrator(int nmeshes)
        : nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
            for (auto &md : collected_) {
                fds_synthetic_turbulence_if_enabled(md->dt, md->t, md->nm);
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

/// Collector for PredFinal sub-graph (CC_IBM only).
///
/// Gathers N kernel results, runs sequential CC_VELOCITY_BC,
/// then emits single BarrierData downstream.
class PredFinalCCCollector
    : public hh::AbstractState<1, MeshData, BarrierData> {
public:
    explicit PredFinalCCCollector(int nmeshes)
        : nmeshes_(nmeshes) {
        results_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        results_.push_back(data);

        if (static_cast<int>(results_.size()) == nmeshes_) {
            std::sort(results_.begin(), results_.end(),
                      [](const auto &a, const auto &b) {
                          return a->nm < b->nm;
                      });

            for (auto &md : results_) {
                fds_cc_velocity_bc(md->t, md->nm, 1);  // applyToEstimated=1 (predictor)
            }

            auto bd = std::make_shared<BarrierData>();
            bd->meshes = std::move(results_);
            results_ = {};
            results_.reserve(nmeshes_);
            this->addResult(bd);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> results_;
};

/// Collector for CorrFinal sub-graph.
///
/// Gathers N kernel results, runs sequential finalization:
///   - CC_VELOCITY_BC (if CC_IBM active — handled internally by fds_cc_velocity_bc)
///   - UPDATE_GLOBAL_OUTPUTS (per-mesh output accumulation — always needed)
/// Then emits single BarrierData downstream (avoids re-collection at graph boundary).
class CorrFinalCollector
    : public hh::AbstractState<1, MeshData, BarrierData> {
public:
    explicit CorrFinalCollector(int nmeshes)
        : nmeshes_(nmeshes) {
        results_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        results_.push_back(data);

        if (static_cast<int>(results_.size()) == nmeshes_) {
            std::sort(results_.begin(), results_.end(),
                      [](const auto &a, const auto &b) {
                          return a->nm < b->nm;
                      });

            for (auto &md : results_) {
                fds_cc_velocity_bc(md->t, md->nm, 0);  // applyToEstimated=0 (corrector)
                fds_update_global_outputs(md->t, md->dt, md->nm);
            }

            auto bd = std::make_shared<BarrierData>();
            bd->meshes = std::move(results_);
            results_ = {};
            results_.reserve(nmeshes_);
            this->addResult(bd);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> results_;
};

#endif // VELOCITY_BC_STATE_H
