// UNUSED — Block decomposition disabled. Kept for reference.
#ifndef WALLBC_BLOCK_SUBGRAPH_H
#define WALLBC_BLOCK_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include <algorithm>
#include <unordered_map>
#include "../data/mesh_data.h"
#include "../data/wallbc_data.h"
#include "../fds_fortran_interface.h"

/// Block kernel task for WallBC.
class WallBCBlockKernelTask
    : public hh::AbstractTask<1, WallBCBlockWork, WallBCBlockWork> {
public:
    explicit WallBCBlockKernelTask(size_t numThreads)
        : hh::AbstractTask<1, WallBCBlockWork, WallBCBlockWork>(
              "WallBCBlockKernel", numThreads) {}

    void execute(std::shared_ptr<WallBCBlockWork> work) override {
        fds_wall_bc_process_cells_block_kernel(
            work->nm, work->t, work->dt, work->dt_bc,
            work->call_ht_1d, work->k1, work->k2);
        this->addResult(work);
    }

    std::shared_ptr<hh::AbstractTask<1, WallBCBlockWork, WallBCBlockWork>>
    copy() override {
        return std::make_shared<WallBCBlockKernelTask>(this->numberThreads());
    }
};

/// Orchestrator task for WallBC block decomposition.
class WallBCBlockOrchestrator
    : public hh::AbstractTask<1, MeshData<>, WallBCBlockWork> {
public:
    WallBCBlockOrchestrator(int nmeshes, int numBlocks)
        : hh::AbstractTask<1, MeshData<>, WallBCBlockWork>("WallBCBlockOrch", 1),
          nmeshes_(nmeshes), numBlocks_(std::max(1, numBlocks)) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData<>> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
            double dt_bc = 0.0;
            int call_ht_1d = 0;

            if (collected_[0]->phase == 1) {
                dt_bc = fds_compute_wall_bc_dt_bc(collected_[0]->t);
                fds_increment_wall_counter();
                call_ht_1d = fds_check_call_ht_1d();
                if (call_ht_1d) {
                    fds_update_bc_clock(collected_[0]->t);
                }
            }

            for (auto &md : collected_) {
                fds_wall_bc_preprocessing_kernel(
                    md->nm, md->t, dt_bc, call_ht_1d);

                int kbar = fds_get_kbar(md->nm);
                int bs = std::max(1, (kbar + numBlocks_ - 1) / numBlocks_);
                int total = (kbar + bs - 1) / bs;

                for (int b = 0; b < total; ++b) {
                    int k1 = b * bs + 1;
                    int k2 = std::min((b + 1) * bs, kbar);
                    this->addResult(std::make_shared<WallBCBlockWork>(
                        md->nm, k1, k2, md->t, md->dt,
                        dt_bc, call_ht_1d, total, md));
                }
            }

            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    int numBlocks_;
    std::vector<std::shared_ptr<MeshData<>>> collected_;
};

/// Collector task for WallBC block decomposition.
class WallBCBlockCollector
    : public hh::AbstractTask<1, WallBCBlockWork, MeshData<>> {
public:
    explicit WallBCBlockCollector(int nmeshes)
        : hh::AbstractTask<1, WallBCBlockWork, MeshData<>>("WallBCBlockColl", 1),
          nmeshes_(nmeshes), nmOffset_(fds_get_lower_mesh_index()) {}

    void execute(std::shared_ptr<WallBCBlockWork> block) override {
        int nm = block->nm;
        auto &entry = entries_[nm];
        if (entry.count == 0) {
            entry.expected = block->totalBlocks;
            entry.meshData = block->originalMeshData;
            entry.dt_bc = block->dt_bc;
            entry.call_ht_1d = block->call_ht_1d;
            entry.isCorrector = (block->originalMeshData->phase == 1);
        }
        entry.count++;

        if (entry.count == entry.expected) {
            completedMeshes_.push_back(entry);
            entries_.erase(nm);

            if (static_cast<int>(completedMeshes_.size()) == nmeshes_) {
                std::sort(completedMeshes_.begin(), completedMeshes_.end(),
                    [this](const Entry &a, const Entry &b) {
                        return a.meshData->nm < b.meshData->nm;
                    });

                for (auto &e : completedMeshes_) {
                    fds_wall_bc_finalize(e.meshData->nm, e.meshData->t,
                                         e.dt_bc, e.call_ht_1d);
                }

                if (completedMeshes_[0].isCorrector) {
                    fds_reset_wall_counter();
                }

                for (auto &e : completedMeshes_) {
                    this->addResult(e.meshData);
                }

                completedMeshes_.clear();
            }
        }
    }

private:
    struct Entry {
        int count = 0;
        int expected = 0;
        double dt_bc = 0.0;
        int call_ht_1d = 0;
        bool isCorrector = false;
        std::shared_ptr<MeshData<>> meshData;
    };
    int nmeshes_;
    int nmOffset_;
    std::unordered_map<int, Entry> entries_;
    std::vector<Entry> completedMeshes_;
};

/// Build the WallBC sub-graph with intra-mesh K-block decomposition.
inline auto buildWallBCBlockSubgraph(int nmeshes, size_t kernelThreads,
                                      int numBlocks) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData<>, MeshData<>>>(
        "WallBCBlock");

    auto orchestratorTask = std::make_shared<WallBCBlockOrchestrator>(nmeshes, numBlocks);
    auto blockKernel = std::make_shared<WallBCBlockKernelTask>(kernelThreads);
    auto collectorTask = std::make_shared<WallBCBlockCollector>(nmeshes);

    subgraph->inputs(orchestratorTask);
    subgraph->edges(orchestratorTask, blockKernel);
    subgraph->edges(blockKernel, collectorTask);
    subgraph->outputs(collectorTask);

    return subgraph;
}

#endif // WALLBC_BLOCK_SUBGRAPH_H
