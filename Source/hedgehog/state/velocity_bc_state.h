#ifndef VELOCITY_BC_STATE_H
#define VELOCITY_BC_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../fds_fortran_interface.h"

/// Collector for PredFinal sub-graph (CC_IBM only).
///
/// Gathers N kernel results, runs sequential CC_VELOCITY_BC,
/// then performs phase transition and emits N MeshData downstream.
class PredFinalCCCollector
    : public hh::AbstractState<1, MeshData, MeshData> {
public:
    explicit PredFinalCCCollector(int nmeshes)
        : nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_[data->nm - nmOffset_] = data;
        ++count_;

        if (count_ == nmeshes_) {
            // CC_VELOCITY_BC finalization
            for (auto &md : collected_) {
                fds_cc_velocity_bc(md->t, md->nm, 1, 1);
            }

            // Phase transition (merged from PhaseTransitionTask)
            double t = collected_[0]->t + collected_[0]->dt;
            double dt = collected_[0]->dt;
            fds_set_predictor(0);
            fds_zero_q_m_dot();
            fds_create_or_remove_obstructions(t, dt);

            count_ = 0;
            for (auto &md : collected_) {
                md->t = t;
                md->phase = 1;
                this->addResult(md);
                md = nullptr;
            }
        }
    }

private:
    int nmeshes_;
    int nmOffset_;
    int count_ = 0;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

/// Collector for PredFinal sub-graph (non-CC_IBM).
///
/// Gathers N kernel results, performs phase transition, emits N MeshData.
/// Merged from former CollectorState + PhaseTransitionTask.
class PredFinalCollector
    : public hh::AbstractState<1, MeshData, MeshData> {
public:
    explicit PredFinalCollector(int nmeshes)
        : nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_[data->nm - nmOffset_] = data;
        ++count_;

        if (count_ == nmeshes_) {
            // Phase transition (merged from PhaseTransitionTask)
            double t = collected_[0]->t + collected_[0]->dt;
            double dt = collected_[0]->dt;
            fds_set_predictor(0);
            fds_zero_q_m_dot();
            fds_create_or_remove_obstructions(t, dt);

            count_ = 0;
            for (auto &md : collected_) {
                md->t = t;
                md->phase = 1;
                this->addResult(md);
                md = nullptr;
            }
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
                fds_cc_velocity_bc(md->t, md->nm, 0, 1);  // applyToEstimated=0 (corrector), DO_IBEDGES=TRUE
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

/// Orchestrator for CorrFinal sub-graph.
///
/// Collects N mesh tokens, runs MeshExchange(6) + CC_END_STEP (merged from
/// former MeshExchange(6b) barrier), then dispatches for parallel kernel.
class CorrFinalOrchestrator
    : public hh::AbstractState<1, MeshData, MeshData> {
public:
    explicit CorrFinalOrchestrator(int nmeshes, bool ccIBM)
        : nmeshes_(nmeshes), ccIBM_(ccIBM) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // Merged from MeshExchange(6b) barrier
            if (ccIBM_) {
                fds_cc_end_step(collected_[0]->t, collected_[0]->dt, 0);
            }
            fds_mesh_exchange(6);

            for (auto &md : collected_) {
                this->addResult(md);
            }
            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    bool ccIBM_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

#endif // VELOCITY_BC_STATE_H
