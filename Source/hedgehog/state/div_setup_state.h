#ifndef DIV_SETUP_STATE_H
#define DIV_SETUP_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Orchestrator state for predictor div setup sub-graph.
/// CC_VELOCITY_BC moved to parallel DivSetupKernelTask (thread-safe).
/// This orchestrator is now a pass-through barrier (kept for graph topology).
class PredDivSetupOrchestrator
    : public hh::AbstractState<1, MeshData, MeshData> {
public:
    explicit PredDivSetupOrchestrator(int nmeshes)
        : hh::AbstractState<1, MeshData, MeshData>(),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
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

/// Orchestrator state for corrector div setup sub-graph.
/// CC_VELOCITY_BC moved to parallel DivSetupKernelTask (thread-safe).
/// This orchestrator is now a pass-through barrier (kept for graph topology).
class CorrDivSetupOrchestrator
    : public hh::AbstractState<1, MeshData, MeshData> {
public:
    explicit CorrDivSetupOrchestrator(int nmeshes)
        : hh::AbstractState<1, MeshData, MeshData>(),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
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

#endif // DIV_SETUP_STATE_H
