#ifndef CORR_PARTICLE_STATE_H
#define CORR_PARTICLE_STATE_H

#include <hedgehog/hedgehog.h>
#include <algorithm>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/corr_particle_data.h"
#include "../fds_fortran_interface.h"

/// Orchestrator for CorrParticle sub-graph (Pattern B).
/// Collects N MeshData tokens, runs sequential PARTICLE_MASS_ENERGY + MOVE_PARTICLES
/// (cross-mesh particle transfer), then dispatches parallel PARTICLE_MOMENTUM_KERNEL.
class CorrParticleOrchestrator : public hh::AbstractState<1, MeshData, CorrParticleWork> {
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
            // Dispatch parallel kernel work
            for (auto &md : collected_) {
                this->addResult(std::make_shared<CorrParticleWork>(md->nm, md->t, md->dt, md));
            }
            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }
};

/// Collector for CorrParticle sub-graph.
class CorrParticleCollector : public hh::AbstractState<1, CorrParticleWork, MeshData> {
    int nmeshes_;
    std::vector<std::shared_ptr<CorrParticleWork>> results_;
public:
    explicit CorrParticleCollector(int nmeshes)
        : nmeshes_(nmeshes) { results_.reserve(nmeshes); }

    void execute(std::shared_ptr<CorrParticleWork> work) override {
        results_.push_back(work);
        if (static_cast<int>(results_.size()) == nmeshes_) {
            std::sort(results_.begin(), results_.end(),
                      [](const auto &a, const auto &b) { return a->nm < b->nm; });
            for (auto &w : results_) {
                this->addResult(w->originalMeshData);
            }
            results_.clear();
            results_.reserve(nmeshes_);
        }
    }
};

#endif // CORR_PARTICLE_STATE_H
