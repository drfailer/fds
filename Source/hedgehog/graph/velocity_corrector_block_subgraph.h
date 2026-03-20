#ifndef VELOCITY_CORRECTOR_BLOCK_SUBGRAPH_H
#define VELOCITY_CORRECTOR_BLOCK_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/mesh_block_data.h"
#include "../state/mesh_block_state.h"
#include "../fds_fortran_interface.h"

/// Merged pre-decompose task for velocity corrector.
/// Combines CC_PROJECT_VELOCITY store + WALL_VELOCITY_NO_GRADH store
/// into a single task to eliminate 1 intermediate queue.
/// Both are no-ops for non-CC_IBM / FFT cases.
class VelCorrPreDecomposeTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit VelCorrPreDecomposeTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "VelCorrPreDecompose", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_cc_project_velocity_kernel(data->nm, data->dt,
                                        /*store=*/1, /*predictor=*/0);
        fds_wall_velocity_no_gradh_kernel(data->nm, data->dt,
                                           /*store=*/1, /*predictor=*/0);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<VelCorrPreDecomposeTask>(this->numberThreads());
    }
};

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

/// Merged post-reassembly task for velocity corrector.
/// Combines CC_PROJECT_VELOCITY fix + WALL_VELOCITY_NO_GRADH fix + CHECK_DIVERGENCE
/// into a single task to eliminate 2 intermediate queues.
/// CC_PROJECT_VELOCITY and WALL_VELOCITY_NO_GRADH are no-ops for non-CC_IBM / FFT cases.
class VelCorrPostReassembleTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit VelCorrPostReassembleTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "VelCorrPostReassemble", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_cc_project_velocity_kernel(data->nm, data->dt,
                                        /*store=*/0, /*predictor=*/0);
        fds_wall_velocity_no_gradh_kernel(data->nm, data->dt,
                                           /*store=*/0, /*predictor=*/0);
        fds_check_divergence_kernel(data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<VelCorrPostReassembleTask>(this->numberThreads());
    }
};

/// Build the velocity corrector sub-graph with block decomposition.
///
/// Pipeline:
///   MeshData -> VelCorrPreDecompose(CCProjectVelStore + WallVelStore)
///            -> Decompose -> VelCorrBlockKernel(parallel) -> Reassemble
///            -> VelCorrPostReassemble(CCProjectVelFix + WallVelFix + CheckDiv)
///            -> MeshData
///
/// @param kernelThreads Number of threads for parallel tasks
/// @param numBlocks Target number of blocks per mesh (default: kernelThreads)
inline auto buildVelocityCorrectorBlockSubgraph(size_t kernelThreads, int numBlocks) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData, MeshData>>("VelocityCorrectorBlock");

    auto preDecompose = std::make_shared<VelCorrPreDecomposeTask>(kernelThreads);
    auto decomposeSM = std::make_shared<hh::StateManager<1, MeshData, MeshBlockData>>(
        std::make_shared<MeshBlockDecomposeState>(numBlocks), "VelCorrDecompose");
    auto blockKernel = std::make_shared<VelocityCorrectorBlockKernelTask>(kernelThreads);
    auto reassembleSM = std::make_shared<hh::StateManager<1, MeshBlockData, MeshData>>(
        std::make_shared<MeshBlockReassembleState>(), "VelCorrReassemble");
    auto postReassemble = std::make_shared<VelCorrPostReassembleTask>(kernelThreads);

    subgraph->inputs(preDecompose);
    subgraph->edges(preDecompose, decomposeSM);
    subgraph->edges(decomposeSM, blockKernel);
    subgraph->edges(blockKernel, reassembleSM);
    subgraph->edges(reassembleSM, postReassemble);
    subgraph->outputs(postReassemble);

    return subgraph;
}

#endif // VELOCITY_CORRECTOR_BLOCK_SUBGRAPH_H
