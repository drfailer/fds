// UNUSED — Block decomposition disabled. Kept for reference.
#ifndef DIVERGENCE_PART2_BLOCK_SUBGRAPH_H
#define DIVERGENCE_PART2_BLOCK_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/mesh_block_data.h"
#include "../state/mesh_block_state.h"
#include "../fds_fortran_interface.h"

/// Block kernel task for DIVERGENCE_PART_2.
/// Processes a K-range sub-block: pressure zone DP, solid zeroing, BC_LOOP, DIV+DDDT.
class DivergencePart2BlockKernelTask
    : public hh::AbstractTask<1, MeshBlockData, MeshBlockData> {
public:
    explicit DivergencePart2BlockKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshBlockData, MeshBlockData>(
              "DivPart2BlockKernel", numThreads) {}

    void execute(std::shared_ptr<MeshBlockData> block) override {
        fds_divergence_part_2_block_kernel(block->nm, block->dt,
                                           block->k1, block->k2);
        this->addResult(block);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshBlockData, MeshBlockData>>
    copy() override {
        return std::make_shared<DivergencePart2BlockKernelTask>(
            this->numberThreads());
    }
};

/// Orchestrator state for DIVERGENCE_PART_2 block decomposition.
/// Collects N MeshData tokens, runs sequential preprocessing per mesh
/// (zone ops, R_PBAR, D_PBAR_DT_P), then decomposes each mesh into K-blocks.
class DivPart2BlockOrchestrator
    : public hh::AbstractState<1, MeshData, MeshBlockData> {
public:
    DivPart2BlockOrchestrator(int nmeshes, int numBlocks)
        : hh::AbstractState<1, MeshData, MeshBlockData>(),
          nmeshes_(nmeshes), numBlocks_(std::max(1, numBlocks)) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // Sequential preprocessing per mesh (zone ops modify global USUM)
            for (auto &md : collected_) {
                fds_divergence_part_2_preprocessing(md->nm, md->dt);
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

/// Build the DIVERGENCE_PART_2 sub-graph with block decomposition.
///
/// Pipeline:
///   MeshData -> Orchestrator(zone ops + decompose) ->
///   DivPart2BlockKernel(parallel) -> Reassemble -> MeshData
///
/// @param nmeshes Number of meshes
/// @param blockThreads Number of threads for parallel block kernel
/// @param numBlocks Target number of blocks per mesh
inline auto buildDivergencePart2BlockSubgraph(int nmeshes, size_t blockThreads,
                                               int numBlocks) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData, MeshData>>("DivPart2Block");

    auto orchestratorSM = std::make_shared<hh::StateManager<1, MeshData, MeshBlockData>>(
        std::make_shared<DivPart2BlockOrchestrator>(nmeshes, numBlocks),
        "DivPart2Orch");
    auto blockKernel = std::make_shared<DivergencePart2BlockKernelTask>(blockThreads);
    auto reassembleSM = std::make_shared<hh::StateManager<1, MeshBlockData, MeshData>>(
        std::make_shared<MeshBlockReassembleState>(), "DivPart2Reassemble");

    subgraph->inputs(orchestratorSM);
    subgraph->edges(orchestratorSM, blockKernel);
    subgraph->edges(blockKernel, reassembleSM);
    subgraph->outputs(reassembleSM);

    return subgraph;
}

#endif // DIVERGENCE_PART2_BLOCK_SUBGRAPH_H
