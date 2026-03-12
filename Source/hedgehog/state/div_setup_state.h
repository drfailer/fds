#ifndef DIV_SETUP_STATE_H
#define DIV_SETUP_STATE_H

#include <algorithm>
#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/div_setup_data.h"
#include "../fds_fortran_interface.h"

/// Orchestrator state for predictor div setup sub-graph.
/// Sequential pre-processing: CC_VELOCITY_BC (if CC_IBM).
/// SET_BAROCLINIC_FALSE, VISCOSITY_BC, AGGLOMERATION moved to parallel kernel task.
class PredDivSetupOrchestrator
    : public hh::AbstractState<1, MeshData, DivSetupWork> {
public:
    explicit PredDivSetupOrchestrator(int nmeshes)
        : hh::AbstractState<1, MeshData, DivSetupWork>(),
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

            // Dispatch parallel kernel work
            for (auto &md : collected_) {
                auto work = std::make_shared<DivSetupWork>(
                    md->nm, md->t, md->dt, 0, md);
                this->addResult(work);
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
    : public hh::AbstractState<1, MeshData, DivSetupWork> {
public:
    explicit CorrDivSetupOrchestrator(int nmeshes)
        : hh::AbstractState<1, MeshData, DivSetupWork>(),
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

            // Dispatch parallel kernel work
            for (auto &md : collected_) {
                auto work = std::make_shared<DivSetupWork>(
                    md->nm, md->t, md->dt, 1, md);
                this->addResult(work);
            }

            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

/// Collector state for div setup sub-graph (shared by predictor and corrector).
class DivSetupCollector
    : public hh::AbstractState<1, DivSetupWork, MeshData> {
public:
    explicit DivSetupCollector(int nmeshes)
        : hh::AbstractState<1, DivSetupWork, MeshData>(),
          nmeshes_(nmeshes) {
        results_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<DivSetupWork> work) override {
        results_.push_back(work);

        if (static_cast<int>(results_.size()) == nmeshes_) {
            std::sort(results_.begin(), results_.end(),
                      [](const auto &a, const auto &b) {
                          return a->nm < b->nm;
                      });

            for (auto &w : results_) {
                this->addResult(w->originalMeshData);
            }

            results_.clear();
            results_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<DivSetupWork>> results_;
};

#endif // DIV_SETUP_STATE_H
