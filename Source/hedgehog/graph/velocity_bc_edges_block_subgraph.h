// UNUSED — Block decomposition disabled. Kept for reference.
#ifndef VELOCITY_BC_EDGES_BLOCK_SUBGRAPH_H
#define VELOCITY_BC_EDGES_BLOCK_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include "../data/mesh_data.h"
#include "../data/barrier_data.h"
#include "../data/velocity_bc_edges_block_data.h"
#include "../fds_fortran_interface.h"

/// Block kernel task for VelocityBC edge processing.
class VelocityBCEdgesBlockKernelTask
    : public hh::AbstractTask<1, VelocityBCEdgesBlockWork, VelocityBCEdgesBlockWork> {
public:
    explicit VelocityBCEdgesBlockKernelTask(size_t numThreads)
        : hh::AbstractTask<1, VelocityBCEdgesBlockWork, VelocityBCEdgesBlockWork>(
              "VelBCEdgesBlockKernel", numThreads) {}

    void execute(std::shared_ptr<VelocityBCEdgesBlockWork> work) override {
        double dragOut = 0.0;
        fds_velocity_bc_process_edges_block_kernel(
            work->nm, work->t, work->applyToEstimated,
            work->k1, work->k2, &dragOut);
        work->dragUvwMax = dragOut;
        this->addResult(work);
    }

    std::shared_ptr<hh::AbstractTask<1, VelocityBCEdgesBlockWork, VelocityBCEdgesBlockWork>>
    copy() override {
        return std::make_shared<VelocityBCEdgesBlockKernelTask>(this->numberThreads());
    }
};

/// Orchestrator task for VelocityBC edges block decomposition.
class VelocityBCEdgesBlockOrchestrator
    : public hh::AbstractTask<1, MeshData<>, VelocityBCEdgesBlockWork> {
public:
    VelocityBCEdgesBlockOrchestrator(int nmeshes, int numBlocks,
                                      int applyToEstimated, bool runSyntheticTurbulence)
        : hh::AbstractTask<1, MeshData<>, VelocityBCEdgesBlockWork>("VelBCEdgesOrch", 1),
          nmeshes_(nmeshes), numBlocks_(std::max(1, numBlocks)),
          applyToEstimated_(applyToEstimated),
          runSyntheticTurbulence_(runSyntheticTurbulence) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData<>> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
            if (runSyntheticTurbulence_) {
                for (auto &md : collected_) {
                    fds_synthetic_turbulence_if_enabled(md->dt, md->t, md->nm);
                }
            }

            for (auto &md : collected_) {
                fds_cc_velocity_cutfaces_ts(md->nm, applyToEstimated_);
                fds_match_velocity_kernel(md->nm, applyToEstimated_);
                fds_velocity_bc_preprocessing(md->nm, md->t, applyToEstimated_);
            }

            for (auto &md : collected_) {
                int kbar = fds_get_kbar(md->nm);
                int bs = std::max(1, (kbar + numBlocks_ - 1) / numBlocks_);
                int total = (kbar + bs - 1) / bs;

                for (int b = 0; b < total; ++b) {
                    int k1 = b * bs + 1;
                    int k2 = std::min((b + 1) * bs, kbar);
                    this->addResult(std::make_shared<VelocityBCEdgesBlockWork>(
                        md->nm, k1, k2, md->t, applyToEstimated_, total, md));
                }
            }

            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    int numBlocks_;
    int applyToEstimated_;
    bool runSyntheticTurbulence_;
    std::vector<std::shared_ptr<MeshData<>>> collected_;
};

/// Collector task for VelocityBC edges block decomposition.
class VelocityBCEdgesBlockCollector
    : public hh::AbstractTask<1, VelocityBCEdgesBlockWork, BarrierData> {
public:
    VelocityBCEdgesBlockCollector(int nmeshes, int applyToEstimated, bool isCorrFinal)
        : hh::AbstractTask<1, VelocityBCEdgesBlockWork, BarrierData>("VelBCEdgesColl", 1),
          nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()),
          applyToEstimated_(applyToEstimated), isCorrFinal_(isCorrFinal) {}

    void execute(std::shared_ptr<VelocityBCEdgesBlockWork> block) override {
        int nm = block->nm;
        auto &entry = entries_[nm];
        if (entry.count == 0) {
            entry.expected = block->totalBlocks;
            entry.meshData = block->originalMeshData;
            entry.dragUvwMax = 0.0;
        }
        entry.dragUvwMax = std::max(entry.dragUvwMax, block->dragUvwMax);
        entry.count++;

        if (entry.count == entry.expected) {
            completedMeshes_.push_back(entry);
            entries_.erase(nm);

            if (static_cast<int>(completedMeshes_.size()) == nmeshes_) {
                std::sort(completedMeshes_.begin(), completedMeshes_.end(),
                    [](const Entry &a, const Entry &b) {
                        return a.meshData->nm < b.meshData->nm;
                    });

                for (auto &e : completedMeshes_) {
                    fds_set_drag_uvwmax(e.meshData->nm, e.dragUvwMax);
                    fds_cc_velocity_bc_ts(e.meshData->t, e.meshData->nm, applyToEstimated_, 1);
                    if (isCorrFinal_) {
                        fds_update_devices_1_ts(e.meshData->t, e.meshData->dt, e.meshData->nm);
                        fds_update_hrr_ts(e.meshData->dt, e.meshData->nm);
                        fds_update_mass_ts(e.meshData->dt, e.meshData->nm);
                        fds_update_fire_spread_outputs_ts(e.meshData->t, e.meshData->dt, e.meshData->nm);
                    }
                }

                if (isCorrFinal_) {
                    fds_reduce_hrr_mass(completedMeshes_[0].meshData->dt);
                }

                auto bd = std::make_shared<BarrierData>();
                bd->meshes.resize(nmeshes_);
                for (auto &e : completedMeshes_) {
                    bd->meshes[e.meshData->nm - nmOffset_] = e.meshData;
                }

                completedMeshes_.clear();
                this->addResult(bd);
            }
        }
    }

private:
    struct Entry {
        int count = 0;
        int expected = 0;
        double dragUvwMax = 0.0;
        std::shared_ptr<MeshData<>> meshData;
    };
    int nmeshes_;
    int nmOffset_;
    int applyToEstimated_;
    bool isCorrFinal_;
    std::unordered_map<int, Entry> entries_;
    std::vector<Entry> completedMeshes_;
};

/// Build the PredFinal sub-graph with block decomposition.
inline auto buildPredFinalBlockSubgraph(int nmeshes, size_t blockThreads, int numBlocks) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData<>, BarrierData>>("PredFinal");

    auto orchTask = std::make_shared<VelocityBCEdgesBlockOrchestrator>(
        nmeshes, numBlocks, /*applyToEstimated=*/1, /*runSyntheticTurbulence=*/true);
    auto blockKernel = std::make_shared<VelocityBCEdgesBlockKernelTask>(blockThreads);
    auto collectorTask = std::make_shared<VelocityBCEdgesBlockCollector>(
        nmeshes, /*applyToEstimated=*/1, /*isCorrFinal=*/false);

    subgraph->inputs(orchTask);
    subgraph->edges(orchTask, blockKernel);
    subgraph->edges(blockKernel, collectorTask);
    subgraph->outputs(collectorTask);

    return subgraph;
}

/// Build the CorrFinal sub-graph with block decomposition.
inline auto buildCorrFinalBlockSubgraph(int nmeshes, size_t blockThreads, int numBlocks) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData<>, BarrierData>>("CorrFinal");

    auto orchTask = std::make_shared<VelocityBCEdgesBlockOrchestrator>(
        nmeshes, numBlocks, /*applyToEstimated=*/0, /*runSyntheticTurbulence=*/false);
    auto blockKernel = std::make_shared<VelocityBCEdgesBlockKernelTask>(blockThreads);
    auto collectorTask = std::make_shared<VelocityBCEdgesBlockCollector>(
        nmeshes, /*applyToEstimated=*/0, /*isCorrFinal=*/true);

    subgraph->inputs(orchTask);
    subgraph->edges(orchTask, blockKernel);
    subgraph->edges(blockKernel, collectorTask);
    subgraph->outputs(collectorTask);

    return subgraph;
}

#endif // VELOCITY_BC_EDGES_BLOCK_SUBGRAPH_H
