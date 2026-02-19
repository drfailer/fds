#ifndef PHASE_TRANSITION_STATE_H
#define PHASE_TRANSITION_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Barrier for the predictor->corrector transition.
/// Collects all mesh tokens from predictor, sets CORRECTOR=TRUE,
/// advances T += DT, creates/removes obstructions, then re-emits
/// tokens with updated time and phase=1 (corrector).
class PhaseTransitionState : public hh::AbstractState<1, MeshData, MeshData> {
public:
    explicit PhaseTransitionState(int nmeshes)
        : hh::AbstractState<1, MeshData, MeshData>(),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nmeshes_) {
            double t = collected_[0]->t;
            double dt = collected_[0]->dt;

            // Transition to corrector phase
            fds_set_predictor(0);  // CORRECTOR=TRUE, PREDICTOR=FALSE

            // Advance time
            t += dt;

            // Zero energy/mass balance arrays
            fds_zero_q_m_dot();

            // Check for obstruction creation/removal
            fds_create_or_remove_obstructions(t, dt);

            // CC_IBM end step (predictor side)
            fds_cc_end_step(t, dt, 0);

            // Re-emit tokens with corrector phase
            for (auto &md : collected_) {
                md->t = t;
                md->phase = 1;
                this->addResult(md);
            }
            collected_.clear();
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

#endif // PHASE_TRANSITION_STATE_H
