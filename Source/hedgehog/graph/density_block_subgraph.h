#ifndef DENSITY_BLOCK_SUBGRAPH_H
#define DENSITY_BLOCK_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/mesh_block_data.h"
#include "../state/mesh_block_state.h"
#include "../fds_fortran_interface.h"

/// Block kernel task for density computation.
/// Processes a K-range sub-block: species mass density, M_DOT_PPP, RHOS/RHO sum.
class DensityBlockKernelTask
    : public hh::AbstractTask<1, MeshBlockData, MeshBlockData> {
public:
    explicit DensityBlockKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshBlockData, MeshBlockData>(
              "DensityBlockKernel", numThreads) {}

    void execute(std::shared_ptr<MeshBlockData> block) override {
        fds_density_block_kernel(block->nm, block->t, block->dt,
                                 block->k1, block->k2);
        this->addResult(block);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshBlockData, MeshBlockData>>
    copy() override {
        return std::make_shared<DensityBlockKernelTask>(this->numberThreads());
    }
};

/// Orchestrator state for density block decomposition.
/// Collects N MeshData tokens, runs sequential preprocessing per mesh
/// (work array setup, wall corrections, settling velocity),
/// then decomposes each mesh into K-blocks for parallel execution.
class DensityBlockOrchestrator
    : public hh::AbstractState<1, MeshData, MeshBlockData> {
public:
    DensityBlockOrchestrator(int nmeshes, int numBlocks)
        : hh::AbstractState<1, MeshData, MeshBlockData>(),
          nmeshes_(nmeshes), numBlocks_(std::max(1, numBlocks)) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // Sequential preprocessing per mesh
            for (auto &md : collected_) {
                fds_density_block_preprocessing(md->nm, md->t, md->dt);
            }

            // Decompose each mesh into K-blocks
            for (auto &md : collected_) {
                int kbar = fds_get_kbar(md->nm);
                int bs = std::max(1, (kbar + numBlocks_ - 1) / numBlocks_);
                int total = (kbar + bs - 1) / bs;

                for (int b = 0; b < total; ++b) {
                    int k1 = b * bs + 1;
                    int k2 = std::min((b + 1) * bs, kbar);
                    this->addResult(std::make_shared<MeshBlockData>(
                        md->nm, k1, k2, md->t, md->dt, md->phase, total, md));
                }
            }

            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    int numBlocks_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

/// Collector state for density block decomposition.
/// Reassembles blocks back into MeshData, runs sequential postprocessing
/// (CHECK_MASS_DENSITY, mass fraction extraction, PBAR, RSUM, TMP).
class DensityBlockCollector
    : public hh::AbstractState<1, MeshBlockData, MeshData> {
public:
    DensityBlockCollector() = default;

    void execute(std::shared_ptr<MeshBlockData> block) override {
        int nm = block->nm;
        auto &entry = entries_[nm];
        if (entry.count == 0) {
            entry.expected = block->totalBlocks;
            entry.meshData = block->originalMeshData;
        }
        entry.count++;
        if (entry.count == entry.expected) {
            auto md = entry.meshData;
            fds_density_block_postprocessing(md->nm, md->t, md->dt);
            this->addResult(md);
            entries_.erase(nm);
        }
    }

private:
    struct Entry {
        int count = 0;
        int expected = 0;
        std::shared_ptr<MeshData> meshData;
    };
    std::unordered_map<int, Entry> entries_;
};

/// Build the density sub-graph with block decomposition.
///
/// Pipeline:
///   MeshData -> Orchestrator(settling vel, work arrays, wall corr, decompose) ->
///   DensityBlockKernel(parallel: species density, M_DOT_PPP, RHOS/RHO sum) ->
///   Collector(STORE_FLUX, CHECK_MASS_DENSITY, ZZ/=RHO, CLIP, PBAR, RSUM, TMP)
///   -> MeshData
///
/// @param nmeshes Number of meshes
/// @param kernelThreads Number of threads for parallel block tasks
/// @param numBlocks Target number of blocks per mesh
inline auto buildDensityBlockSubgraph(int nmeshes, size_t kernelThreads,
                                       int numBlocks) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData, MeshData>>("DensityBlock");

    auto orchestratorSM = std::make_shared<hh::StateManager<1, MeshData, MeshBlockData>>(
        std::make_shared<DensityBlockOrchestrator>(nmeshes, numBlocks),
        "DensityOrch");
    auto blockKernel = std::make_shared<DensityBlockKernelTask>(kernelThreads);
    auto collectorSM = std::make_shared<hh::StateManager<1, MeshBlockData, MeshData>>(
        std::make_shared<DensityBlockCollector>(), "DensityCollector");

    subgraph->inputs(orchestratorSM);
    subgraph->edges(orchestratorSM, blockKernel);
    subgraph->edges(blockKernel, collectorSM);
    subgraph->outputs(collectorSM);

    return subgraph;
}

#endif // DENSITY_BLOCK_SUBGRAPH_H
