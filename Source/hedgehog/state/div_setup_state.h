#ifndef DIV_SETUP_STATE_H
#define DIV_SETUP_STATE_H

#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Orchestrator state for predictor div setup sub-graph.
/// Sequential pre-processing: CC_VELOCITY_BC (if CC_IBM).
/// SET_BAROCLINIC_FALSE, VISCOSITY_BC, AGGLOMERATION moved to parallel kernel task.
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
            // Sequential pre-processing (only CC_IBM remains)
            for (auto &md : collected_) {
                fds_cc_velocity_bc(md->t, md->nm, 0);  // CC_IBM: sequential (pending conversion)
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

/// Orchestrator state for corrector div setup sub-graph.
/// Sequential pre-processing: CC_VELOCITY_BC (if CC_IBM).
/// SET_BAROCLINIC_FALSE, VISCOSITY_BC, AGGLOMERATION moved to parallel kernel task.
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
            // Sequential pre-processing (only CC_IBM remains)
            for (auto &md : collected_) {
                fds_cc_velocity_bc(md->t, md->nm, 1);  // CC_IBM: sequential (pending conversion)
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

#endif // DIV_SETUP_STATE_H
