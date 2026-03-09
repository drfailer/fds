#ifndef CORR_DIV_PART1_STATE_H
#define CORR_DIV_PART1_STATE_H

#include <algorithm>
#include <hedgehog/hedgehog.h>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/corr_div_part1_data.h"
#include "../fds_fortran_interface.h"

/// Orchestrator state for corrector divergence part 1 sub-graph.
/// Runs sequential COMBUSTION_BC (cross-mesh dependency: reads OMESH%Q)
/// on each mesh before dispatching parallel DIVERGENCE_PART_1_KERNEL work.
class CorrDivPart1Orchestrator
    : public hh::AbstractState<1, MeshData, CorrDivPart1Work> {
public:
    explicit CorrDivPart1Orchestrator(int nmeshes)
        : hh::AbstractState<1, MeshData, CorrDivPart1Work>(),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // Sequential pre-processing: COMBUSTION_BC reads OMESH(NOM)%Q
            for (auto &md : collected_) {
                fds_combustion_bc(md->nm);
            }

            // Dispatch parallel kernel work
            for (auto &md : collected_) {
                auto work = std::make_shared<CorrDivPart1Work>(
                    md->nm, md->t, md->dt, md);
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

/// Collector state for corrector divergence part 1 sub-graph.
class CorrDivPart1Collector
    : public hh::AbstractState<1, CorrDivPart1Work, MeshData> {
public:
    explicit CorrDivPart1Collector(int nmeshes)
        : hh::AbstractState<1, CorrDivPart1Work, MeshData>(),
          nmeshes_(nmeshes) {
        results_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<CorrDivPart1Work> work) override {
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
    std::vector<std::shared_ptr<CorrDivPart1Work>> results_;
};

#endif // CORR_DIV_PART1_STATE_H
