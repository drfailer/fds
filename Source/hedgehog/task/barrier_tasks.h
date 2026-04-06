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
// Tasks that receive BarrierData from upstream collectors and scatter MeshData
// downstream.  Most barrier patterns use BarrierState (state/barrier_state.h)
// which merges collector + barrier into a single MeshData→MeshData state node.
// ---------------------------------------------------------------------------

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
        }
        this->batchAddResult(data->meshes);
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

/// CorrFinal dump task — reduces HRR/MASS, checks dump schedule,
/// emits MeshData (per-mesh dump) + BarrierData (global dump + timestep loop).
///
/// Replaces CorrFinalCollector state (which lacked canTerminate).
/// Upstream: CollectorState collects N MeshData → 1 BarrierData.
class CorrFinalDumpTask
    : public hh::AbstractTask<1, BarrierData, MeshData, BarrierData> {
public:
    CorrFinalDumpTask()
        : hh::AbstractTask<1, BarrierData, MeshData, BarrierData>(
              "CorrFinalDump", 1) {}

    void execute(std::shared_ptr<BarrierData> data) override {
        auto t0 = std::chrono::steady_clock::now();

        // Reduce per-mesh accumulators into globals
        fds_reduce_hrr_mass(data->dt());

        // Check if any mesh needs dump I/O this timestep
        bool anyDump = false;
        for (auto &md : data->meshes) {
            bool dump = false;
            fds_check_dump_schedule(md->t, md->nm, &dump);
            if (dump) { anyDump = true; break; }
        }

        // Build BarrierData (always emitted for DumpGlobalTask)
        auto bd = std::make_shared<BarrierData>();
        bd->meshes = data->meshes;  // copy shared_ptrs (TimestepState needs them)

        if (anyDump) {
            bd->skipMeshDump = false;
            // Emit MeshData tokens for per-mesh dump I/O
            this->batchAddResult(data->meshes);
        } else {
            bd->skipMeshDump = true;
        }

        // Always emit BarrierData for DumpGlobalTask
        this->addResult(bd);

        auto t1 = std::chrono::steady_clock::now();
        totalTime_ += std::chrono::duration<double>(t1 - t0).count();
        ++invocations_;
    }

    [[nodiscard]] std::string extraPrintingInformation() const override {
        std::ostringstream oss;
        oss << "REDUCE_HRR_MASS\\n"
            << "CHECK_DUMP_SCHEDULE\\n"
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
