#ifndef BARRIER_TASKS_H
#define BARRIER_TASKS_H

#include <hedgehog/hedgehog.h>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <string>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../fds_fortran_interface.h"

// ---------------------------------------------------------------------------
// Barrier computation tasks.
//
// Remaining tasks that receive BarrierData from upstream collectors or
// sub-graphs and scatter MeshData downstream.  Most former barrier tasks
// have been replaced by BarrierState (state/barrier_state.h) which merges
// the collector + barrier into a single state node.
// ---------------------------------------------------------------------------

/// MESH_EXCHANGE barrier task.
///
/// Used for exchanges where the upstream already emits BarrierData
/// (e.g. MeshExchange(2) after CorrRadiation or Fork2 join).
///
/// Optional pre/post-exchange operations:
/// - ccDensity: run CC_DENSITY(T,DT) before the exchange
/// - ccEndStep: run CC_END_STEP(T,DT) before the exchange
/// - initDiv: run INITIALIZE_DIVERGENCE_INTEGRALS after the exchange
class MeshExchangeTask : public hh::AbstractTask<1, BarrierData, MeshData> {
public:
    explicit MeshExchangeTask(int code, bool ccDensity = false,
                              bool ccEndStep = false, bool initDiv = false)
        : hh::AbstractTask<1, BarrierData, MeshData>(
              "MeshExchange(" + std::to_string(code) + ")", 1),
          code_(code), ccDensity_(ccDensity), ccEndStep_(ccEndStep),
          initDiv_(initDiv) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        auto t0 = std::chrono::steady_clock::now();
        if (ccDensity_) { fds_cc_density(data->t(), data->dt()); }
        if (ccEndStep_) { fds_cc_end_step(data->t(), data->dt(), 0); }
        if (code_ != 2 || fds_exchange_radiation()) { fds_mesh_exchange(code_); }
        if (code_ == 1) { fds_exchange_inserted_particles(); }
        if (initDiv_) { fds_initialize_divergence_integrals(); }
        auto t1 = std::chrono::steady_clock::now();
        totalTime_ += std::chrono::duration<double>(t1 - t0).count();
        ++invocations_;
        for (auto &md : data->meshes) { this->addResult(md); }
    }

    std::string extraPrintingInformation() const override {
        std::ostringstream oss;
        if (ccDensity_) oss << "CC_DENSITY\\n";
        if (ccEndStep_) oss << "CC_END_STEP\\n";
        oss << "MESH_EXCHANGE(" << code_ << ")\\n";
        if (code_ == 1) oss << "EXCHANGE_INSERTED_PARTICLES\\n";
        if (initDiv_) oss << "INITIALIZE_DIVERGENCE_INTEGRALS\\n";
        oss << std::fixed << std::setprecision(3) << totalTime_ << "s"
            << " / " << invocations_ << " calls";
        if (invocations_ > 0)
            oss << " / avg " << std::setprecision(3)
                << (totalTime_ * 1000.0 / invocations_) << "ms";
        return oss.str();
    }

private:
    int code_;
    bool ccDensity_;
    bool ccEndStep_;
    bool initDiv_;
    double totalTime_ = 0.0;
    int invocations_ = 0;
};

/// Merged MeshExchange(2) + DivExchange barrier task (Opt 3).
///
/// Takes BarrierData from Fork2 join, runs MESH_EXCHANGE(2) + QR_ADD per mesh +
/// exchange divergence info + RTE source correction + global matrix reassign,
/// then scatters MeshData downstream.  Replaces two separate nodes in the
/// non-CC_IBM corrector path.
class CorrMeshExch2DivExchangeTask : public hh::AbstractTask<1, BarrierData, MeshData> {
public:
    CorrMeshExch2DivExchangeTask()
        : hh::AbstractTask<1, BarrierData, MeshData>("MeshExch2+DivExch", 1) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        auto t0 = std::chrono::steady_clock::now();
        if (fds_exchange_radiation()) { fds_mesh_exchange(2); }
        for (auto &md : data->meshes) {
            fds_divergence_part_1_add_qr_b(md->nm);
        }
        fds_exchange_divergence_info();
        fds_rte_source_correction();
        fds_global_matrix_reassign(0);
        auto t1 = std::chrono::steady_clock::now();
        totalTime_ += std::chrono::duration<double>(t1 - t0).count();
        ++invocations_;
        for (auto &md : data->meshes) { this->addResult(md); }
    }

    std::string extraPrintingInformation() const override {
        std::ostringstream oss;
        oss << "MESH_EXCHANGE(2)\\n"
            << "QR_ADD\\n"
            << "EXCH_DIV_INFO\\n"
            << "RTE_SOURCE_CORR\\n"
            << "GLOBAL_MATRIX_REASSIGN\\n"
            << std::fixed << std::setprecision(3) << totalTime_ << "s"
            << " / " << invocations_ << " calls";
        if (invocations_ > 0)
            oss << " / avg " << std::setprecision(3)
                << (totalTime_ * 1000.0 / invocations_) << "ms";
        return oss.str();
    }

private:
    double totalTime_ = 0.0;
    int invocations_ = 0;
};

/// Phase transition task — sets CORRECTOR=TRUE, advances T, zeros arrays,
/// handles obstructions.
class PhaseTransitionTask : public hh::AbstractTask<1, BarrierData, MeshData> {
public:
    PhaseTransitionTask()
        : hh::AbstractTask<1, BarrierData, MeshData>("PhaseTransition", 1) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        auto t0 = std::chrono::steady_clock::now();
        double t = data->t();
        double dt = data->dt();

        fds_set_predictor(0);  // CORRECTOR=TRUE, PREDICTOR=FALSE
        t += dt;
        fds_zero_q_m_dot();
        fds_create_or_remove_obstructions(t, dt);

        auto t1 = std::chrono::steady_clock::now();
        totalTime_ += std::chrono::duration<double>(t1 - t0).count();
        ++invocations_;

        for (auto &md : data->meshes) {
            md->t = t;
            md->phase = 1;  // corrector
            this->addResult(md);
        }
    }

    std::string extraPrintingInformation() const override {
        std::ostringstream oss;
        oss << "SET_PREDICTOR(0)\\n"
            << "ZERO_Q_M_DOT\\n"
            << "CREATE_OR_REMOVE_OBSTRUCTIONS\\n"
            << std::fixed << std::setprecision(3) << totalTime_ << "s"
            << " / " << invocations_ << " calls";
        if (invocations_ > 0)
            oss << " / avg " << std::setprecision(3)
                << (totalTime_ * 1000.0 / invocations_) << "ms";
        return oss.str();
    }

private:
    double totalTime_ = 0.0;
    int invocations_ = 0;
};

#endif // BARRIER_TASKS_H
