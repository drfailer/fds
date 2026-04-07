#ifndef CORR_PARTICLE_STATE_H
#define CORR_PARTICLE_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Orchestrator task for CorrParticle sub-graph.
/// Collects N MeshData tokens, runs sequential PARTICLE_MASS_ENERGY + MOVE_PARTICLES
/// (cross-mesh particle transfer), then dispatches MeshData for parallel PARTICLE_MOMENTUM_KERNEL.
///
/// Runs on a single thread.
class CorrParticleOrchestrator : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit CorrParticleOrchestrator(int nmeshes)
        : hh::AbstractTask<1, MeshData, MeshData>("CorrParticleOrch", 1),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // Sequential pre-processing: particle routines (cross-mesh transfer)
            for (auto &md : collected_) {
                fds_particle_mass_energy(md->t, md->dt, md->nm);
                fds_move_particles(md->t, md->dt, md->nm);
            }
            // Dispatch for parallel kernel execution
            this->batchAddResult(collected_);
            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

#endif // CORR_PARTICLE_STATE_H
