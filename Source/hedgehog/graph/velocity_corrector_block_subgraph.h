#ifndef VELOCITY_CORRECTOR_BLOCK_SUBGRAPH_H
#define VELOCITY_CORRECTOR_BLOCK_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/mesh_block_data.h"
#include "../state/mesh_block_state.h"
#include "../fds_fortran_interface.h"

/// Block kernel task for velocity corrector.
/// Processes a K-range sub-block of a single mesh.
class VelocityCorrectorBlockKernelTask
    : public hh::AbstractTask<1, MeshBlockData, MeshBlockData> {
public:
    explicit VelocityCorrectorBlockKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshBlockData, MeshBlockData>(
              "VelCorrBlockKernel", numThreads) {}

    void execute(std::shared_ptr<MeshBlockData> block) override {
        fds_velocity_corrector_block_kernel(block->nm, block->dt,
                                            block->k1, block->k2);
        this->addResult(block);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshBlockData, MeshBlockData>>
    copy() override {
        return std::make_shared<VelocityCorrectorBlockKernelTask>(
            this->numberThreads());
    }
};

/// Mesh-level check divergence task (runs after block reassembly).
class CheckDivergenceKernelTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit CheckDivergenceKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "CheckDivKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_check_divergence_kernel(data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<CheckDivergenceKernelTask>(
            this->numberThreads());
    }
};

/// Build the velocity corrector sub-graph with block decomposition.
///
/// Pipeline:
///   MeshData -> Decompose -> VelCorrBlockKernel(parallel) -> Reassemble -> CheckDivKernel(parallel) -> MeshData
///
/// The velocity corrector kernel (pure I,J,K loops) is block-decomposed along K.
/// The check divergence kernel (reduction to mesh scalars) runs at mesh level after reassembly.
///
/// @param kernelThreads Number of threads for parallel tasks
/// @param numBlocks Target number of blocks per mesh (default: kernelThreads)
inline auto buildVelocityCorrectorBlockSubgraph(size_t kernelThreads, int numBlocks) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData, MeshData>>("VelocityCorrectorBlock");

    auto decomposeSM = std::make_shared<hh::StateManager<1, MeshData, MeshBlockData>>(
        std::make_shared<MeshBlockDecomposeState>(numBlocks), "VelCorrDecompose");
    auto blockKernel = std::make_shared<VelocityCorrectorBlockKernelTask>(kernelThreads);
    auto reassembleSM = std::make_shared<hh::StateManager<1, MeshBlockData, MeshData>>(
        std::make_shared<MeshBlockReassembleState>(), "VelCorrReassemble");
    auto checkDivTask = std::make_shared<CheckDivergenceKernelTask>(kernelThreads);

    subgraph->inputs(decomposeSM);
    subgraph->edges(decomposeSM, blockKernel);
    subgraph->edges(blockKernel, reassembleSM);
    subgraph->edges(reassembleSM, checkDivTask);
    subgraph->outputs(checkDivTask);

    return subgraph;
}

#endif // VELOCITY_CORRECTOR_BLOCK_SUBGRAPH_H
