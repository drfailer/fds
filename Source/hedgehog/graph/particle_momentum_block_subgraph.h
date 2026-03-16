#ifndef PARTICLE_MOMENTUM_BLOCK_SUBGRAPH_H
#define PARTICLE_MOMENTUM_BLOCK_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/mesh_block_data.h"
#include "../state/mesh_block_state.h"
#include "../fds_fortran_interface.h"

/// Block kernel task for particle momentum transfer.
/// Processes a K-range sub-block of a single mesh.
class ParticleMomentumBlockKernelTask
    : public hh::AbstractTask<1, MeshBlockData, MeshBlockData> {
public:
    explicit ParticleMomentumBlockKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshBlockData, MeshBlockData>(
              "PartMomBlockKernel", numThreads) {}

    void execute(std::shared_ptr<MeshBlockData> block) override {
        fds_particle_momentum_block_kernel(block->nm, block->dt,
                                            block->k1, block->k2);
        this->addResult(block);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshBlockData, MeshBlockData>>
    copy() override {
        return std::make_shared<ParticleMomentumBlockKernelTask>(
            this->numberThreads());
    }
};

/// Build the particle momentum sub-graph with block decomposition.
///
/// Pipeline:
///   MeshData -> Decompose -> PartMomBlockKernel(parallel) -> Reassemble -> MeshData
///
/// @param kernelThreads Number of threads for parallel tasks
/// @param numBlocks Target number of blocks per mesh
inline auto buildParticleMomentumBlockSubgraph(size_t kernelThreads, int numBlocks) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData, MeshData>>("ParticleMomentumBlock");

    auto decomposeSM = std::make_shared<hh::StateManager<1, MeshData, MeshBlockData>>(
        std::make_shared<MeshBlockDecomposeState>(numBlocks), "PartMomDecompose");
    auto blockKernel = std::make_shared<ParticleMomentumBlockKernelTask>(kernelThreads);
    auto reassembleSM = std::make_shared<hh::StateManager<1, MeshBlockData, MeshData>>(
        std::make_shared<MeshBlockReassembleState>(), "PartMomReassemble");

    subgraph->inputs(decomposeSM);
    subgraph->edges(decomposeSM, blockKernel);
    subgraph->edges(blockKernel, reassembleSM);
    subgraph->outputs(reassembleSM);

    return subgraph;
}

#endif // PARTICLE_MOMENTUM_BLOCK_SUBGRAPH_H
