#ifndef PRED_STEP1_STATE_H
#define PRED_STEP1_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Orchestrator task for PredStep1 sub-graph.
/// Collects N MeshData tokens, runs sequential INSERT_ALL_PARTICLES for each mesh,
/// then dispatches MeshData for parallel COMPUTE_VISCOSITY + MASS_FINITE_DIFFERENCES kernels.
///
/// Runs on a single thread.
class PredStep1Orchestrator : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit PredStep1Orchestrator(int nmeshes)
        : hh::AbstractTask<1, MeshData, MeshData>("PredStep1Orch", 1),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);
        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // Sequential pre-processing: INSERT_ALL_PARTICLES (cross-mesh, global state)
            for (auto &md : collected_) {
                fds_insert_particles(md->t, md->nm);
            }
            // Dispatch for parallel kernel execution
            for (auto &md : collected_) {
                this->addResult(md);
            }
            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

#endif // PRED_STEP1_STATE_H
