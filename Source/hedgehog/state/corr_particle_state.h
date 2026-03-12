#ifndef CORR_PARTICLE_STATE_H
#define CORR_PARTICLE_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Orchestrator for CorrParticle sub-graph (Pattern B).
/// Collects N MeshData tokens, runs sequential PARTICLE_MASS_ENERGY + MOVE_PARTICLES
/// (cross-mesh particle transfer), then dispatches MeshData for parallel PARTICLE_MOMENTUM_KERNEL.
class CorrParticleOrchestrator : public hh::AbstractState<1, MeshData, MeshData> {
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> collected_;
public:
    explicit CorrParticleOrchestrator(int nmeshes)
        : nmeshes_(nmeshes) { collected_.reserve(nmeshes); }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // Sequential pre-processing: particle routines (cross-mesh transfer)
            for (auto &md : collected_) {
                fds_particle_mass_energy(md->t, md->dt, md->nm);
                fds_move_particles(md->t, md->dt, md->nm);
            }
            // Dispatch for parallel kernel execution
            for (auto &md : collected_) {
                this->addResult(md);
            }
            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }
};

#endif // CORR_PARTICLE_STATE_H
