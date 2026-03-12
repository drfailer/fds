#ifndef VELOCITY_BC_STATE_H
#define VELOCITY_BC_STATE_H

#include <algorithm>
#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/velocity_bc_data.h"
#include "../fds_fortran_interface.h"

/// Orchestrator for PredFinal sub-graph (Pattern B - Sequential Pre-Processing).
///
/// Collects N mesh tokens, runs sequential preprocessing for each mesh:
///   - SYNTHETIC_TURBULENCE_IF_ENABLED (SEM inflow — uses RANDOM_NUMBER, kept sequential)
/// Then dispatches N VelocityBCWork tokens for parallel processing:
///   - MATCH_VELOCITY_KERNEL (thread-safe cross-mesh interpolation via M%)
///   - VELOCITY_BC_PREPROCESSING (thread-safe OMESH reads via M%)
///   - VELOCITY_BC_PROCESS_EDGES_KERNEL (edge boundary conditions)
class PredFinalOrchestrator
    : public hh::AbstractState<1, MeshData, VelocityBCWork> {
public:
    explicit PredFinalOrchestrator(int nmeshes)
        : hh::AbstractState<1, MeshData, VelocityBCWork>(),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // Sequential: SYNTHETIC_TURBULENCE only (RANDOM_NUMBER not thread-safe)
            for (auto &md : collected_) {
                fds_synthetic_turbulence_if_enabled(md->dt, md->t, md->nm);
            }

            // Dispatch parallel work (match_velocity + preprocessing + edges)
            for (auto &md : collected_) {
                auto work = std::make_shared<VelocityBCWork>(
                    md->nm, md->t, md->dt, /*estimated=*/1, md);
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

/// Collector for PredFinal sub-graph (Pattern B - Sequential Post-Processing).
///
/// Gathers N kernel results, runs sequential finalization:
///   - CC_VELOCITY_BC (if CC_IBM is active)
/// Then emits N MeshData tokens downstream.
class PredFinalCollector
    : public hh::AbstractState<1, VelocityBCWork, MeshData> {
public:
    explicit PredFinalCollector(int nmeshes)
        : hh::AbstractState<1, VelocityBCWork, MeshData>(),
          nmeshes_(nmeshes) {
        results_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<VelocityBCWork> work) override {
        results_.push_back(work);

        if (static_cast<int>(results_.size()) == nmeshes_) {
            std::sort(results_.begin(), results_.end(),
                      [](const auto &a, const auto &b) {
                          return a->nm < b->nm;
                      });

            // Sequential finalization
            for (auto &w : results_) {
                fds_cc_velocity_bc(w->t, w->nm, w->applyToEstimated);
            }

            for (auto &w : results_) {
                this->addResult(w->originalMeshData);
            }

            results_.clear();
            results_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<VelocityBCWork>> results_;
};

/// Orchestrator for CorrFinal sub-graph (Pattern B - Sequential Pre-Processing).
///
/// Collects N mesh tokens, then dispatches N VelocityBCWork tokens for parallel processing:
///   - MATCH_VELOCITY_KERNEL (thread-safe cross-mesh interpolation via M%)
///   - VELOCITY_BC_PREPROCESSING (thread-safe OMESH reads via M%)
///   - VELOCITY_BC_PROCESS_EDGES_KERNEL (edge boundary conditions)
/// No sequential preprocessing remains in the corrector orchestrator.
class CorrFinalOrchestrator
    : public hh::AbstractState<1, MeshData, VelocityBCWork> {
public:
    explicit CorrFinalOrchestrator(int nmeshes)
        : hh::AbstractState<1, MeshData, VelocityBCWork>(),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // Dispatch parallel work (match_velocity + preprocessing + edges)
            for (auto &md : collected_) {
                auto work = std::make_shared<VelocityBCWork>(
                    md->nm, md->t, md->dt, /*estimated=*/0, md);
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

/// Collector for CorrFinal sub-graph (Pattern B - Sequential Post-Processing).
///
/// Gathers N kernel results, runs sequential finalization:
///   - CC_VELOCITY_BC (if CC_IBM is active)
///   - UPDATE_GLOBAL_OUTPUTS (per-mesh output accumulation)
/// Then emits N MeshData tokens downstream.
class CorrFinalCollector
    : public hh::AbstractState<1, VelocityBCWork, MeshData> {
public:
    explicit CorrFinalCollector(int nmeshes)
        : hh::AbstractState<1, VelocityBCWork, MeshData>(),
          nmeshes_(nmeshes) {
        results_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<VelocityBCWork> work) override {
        results_.push_back(work);

        if (static_cast<int>(results_.size()) == nmeshes_) {
            std::sort(results_.begin(), results_.end(),
                      [](const auto &a, const auto &b) {
                          return a->nm < b->nm;
                      });

            // Sequential finalization
            for (auto &w : results_) {
                fds_cc_velocity_bc(w->t, w->nm, w->applyToEstimated);
                fds_update_global_outputs(w->t, w->originalMeshData->dt, w->nm);
            }

            for (auto &w : results_) {
                this->addResult(w->originalMeshData);
            }

            results_.clear();
            results_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<VelocityBCWork>> results_;
};

#endif // VELOCITY_BC_STATE_H
