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
// Tasks that receive BarrierData from upstream collectors and scatter MeshData<>
// downstream.  Most barrier patterns use BarrierState (state/barrier_state.h)
// which merges collector + barrier into a single MeshData<>→MeshData<> state node.
// ---------------------------------------------------------------------------

/// Merged collector + phase transition: collects N MeshData<>, runs phase
/// transition (CORRECTOR=TRUE, advance T, zero arrays, obstructions), then
/// re-emits N MeshData<> with updated t and phase=1.
/// Replaces separate CollectorTask + PhaseTransitionTask.
class PhaseTransitionTask : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit PhaseTransitionTask(int nmeshes)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>("PhaseTransition", 1),
          nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()) {
        collected_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshData<>> data) override {
        collected_[data->nm - nmOffset_] = data;
        if (++count_ == nmeshes_) {
            auto t0 = std::chrono::steady_clock::now();
            double t = collected_[0]->t;
            double dt = collected_[0]->dt;

            fds_set_predictor(0);  // CORRECTOR=TRUE, PREDICTOR=FALSE
            t += dt;
            fds_zero_q_m_dot();
            fds_create_or_remove_obstructions(t, dt);

            // Increment WALL_COUNTER once per corrector step, carry per-mesh
            fds_increment_wall_counter();
            int wc = fds_get_wall_counter();

            auto t1 = std::chrono::steady_clock::now();
            totalTime_ += std::chrono::duration<double>(t1 - t0).count();
            ++invocations_;

            for (auto &md : collected_) {
                md->t = t;
                md->phase = 1;  // corrector
                md->wall_counter = wc;
            }
            this->batchAddResult(collected_);
            for (auto &md : collected_) { md = nullptr; }
            count_ = 0;
        }
    }

    [[nodiscard]] std::string extraPrintingInformation() const override {
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
    int nmeshes_, nmOffset_, count_ = 0;
    double totalTime_ = 0.0;
    int invocations_ = 0;
    std::vector<std::shared_ptr<MeshData<>>> collected_;
};

#endif // BARRIER_TASKS_H
