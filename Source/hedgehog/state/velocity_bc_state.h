#ifndef VELOCITY_BC_STATE_H
#define VELOCITY_BC_STATE_H

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
        : nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_[data->nm - nmOffset_] = data;
        ++count_;

        if (count_ == nmeshes_) {
            for (auto &md : collected_) {
                fds_cc_velocity_bc(md->t, md->nm, 1);  // applyToEstimated=1 (predictor)
            }

            auto bd = std::make_shared<BarrierData>();
            bd->meshes = std::move(collected_);
            collected_.resize(nmeshes_, nullptr);
            count_ = 0;
            this->addResult(bd);
        }
    }

private:
    int nmeshes_;
    int nmOffset_;
    int count_ = 0;
    std::vector<std::shared_ptr<MeshData>> collected_;
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
        : nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_[data->nm - nmOffset_] = data;
        ++count_;

        if (count_ == nmeshes_) {
            for (auto &md : collected_) {
                fds_cc_velocity_bc(md->t, md->nm, 0);  // applyToEstimated=0 (corrector)
                fds_update_global_outputs(md->t, md->dt, md->nm);
            }

            auto bd = std::make_shared<BarrierData>();
            bd->meshes = std::move(collected_);
            collected_.resize(nmeshes_, nullptr);
            count_ = 0;
            this->addResult(bd);
        }
    }

private:
    int nmeshes_;
    int nmOffset_;
    int count_ = 0;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

#endif // VELOCITY_BC_STATE_H
