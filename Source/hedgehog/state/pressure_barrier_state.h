#ifndef PRESSURE_BARRIER_STATE_H
#define PRESSURE_BARRIER_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Barrier that wraps the entire pressure iteration scheme.
/// Collects all mesh tokens, calls the opaque Fortran pressure solver,
/// then re-emits all tokens.
class PressureBarrierState : public hh::AbstractState<1, MeshData, MeshData> {
public:
    explicit PressureBarrierState(int nmeshes)
        : hh::AbstractState<1, MeshData, MeshData>(),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // Use time from first token (all should have the same T, DT)
            double t = collected_[0]->t;
            double dt = collected_[0]->dt;
            fds_pressure_iteration(t, dt);
            for (auto &md : collected_) {
                this->addResult(md);
            }
            collected_.clear();
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

#endif // PRESSURE_BARRIER_STATE_H
