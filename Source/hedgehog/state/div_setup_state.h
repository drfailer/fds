#ifndef DIV_SETUP_STATE_H
#define DIV_SETUP_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Orchestrator task for predictor div setup sub-graph.
/// CC_VELOCITY_BC moved to parallel DivSetupKernelTask (thread-safe).
/// This orchestrator is now a pass-through barrier (kept for graph topology).
///
/// Runs on a single thread.
class PredDivSetupOrchestrator
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit PredDivSetupOrchestrator(int nmeshes)
        : hh::AbstractTask<1, MeshData, MeshData>("PredDivSetupOrch", 1),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
            this->batchAddResult(collected_);
            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

/// Orchestrator task for corrector div setup sub-graph.
/// CC_VELOCITY_BC moved to parallel DivSetupKernelTask (thread-safe).
/// This orchestrator is now a pass-through barrier (kept for graph topology).
///
/// Runs on a single thread.
class CorrDivSetupOrchestrator
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit CorrDivSetupOrchestrator(int nmeshes)
        : hh::AbstractTask<1, MeshData, MeshData>("Fork1VFluxOrch", 1),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
            this->batchAddResult(collected_);
            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

#endif // DIV_SETUP_STATE_H
