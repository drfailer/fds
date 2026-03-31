#ifndef VELOCITY_BC_STATE_H
#define VELOCITY_BC_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../fds_fortran_interface.h"

/// Collector for PredFinal sub-graph.
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

/// Collector for CorrFinal sub-graph (merged with PreDumpScatter).
///
/// Gathers N kernel results from VelocityBCEdgesTask (which already ran
/// UPDATE_HRR_TS, UPDATE_MASS_TS, UPDATE_FIRE_SPREAD_OUTPUTS_TS in parallel),
/// then performs:
///   1. REDUCE_HRR_MASS — sum per-mesh Q_DOT_MESH/M_DOT_MESH into globals
///   2. Check dump schedule for all meshes
///   3. Emit N MeshData → DumpMeshOutputsTask (only if any mesh needs dump)
///   4. Emit 1 BarrierData → DumpGlobalTask (always, carries meshes + skipMeshDump flag)
class CorrFinalCollector
    : public hh::AbstractState<1, MeshData, MeshData, BarrierData> {
public:
    explicit CorrFinalCollector(int nmeshes)
        : nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_[data->nm - nmOffset_] = data;
        ++count_;

        if (count_ == nmeshes_) {
            // Reduce per-mesh accumulators into globals
            fds_reduce_hrr_mass(collected_[0]->dt);

            // Check if any mesh needs dump I/O this timestep
            bool anyDump = false;
            for (auto &md : collected_) {
                bool dump = false;
                fds_check_dump_schedule(md->t, md->nm, &dump);
                if (dump) { anyDump = true; break; }
            }

            // Build BarrierData (always emitted for DumpGlobalTask)
            auto bd = std::make_shared<BarrierData>();
            bd->meshes = collected_;  // copy shared_ptrs (TimestepState needs them)

            if (anyDump) {
                bd->skipMeshDump = false;
                // Emit MeshData tokens for per-mesh dump I/O
                for (auto &md : collected_) {
                    this->addResult(md);
                }
            } else {
                bd->skipMeshDump = true;
            }

            // Always emit BarrierData for DumpGlobalTask
            this->addResult(bd);

            // Reset for next timestep
            count_ = 0;
            for (auto &md : collected_) { md = nullptr; }
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
