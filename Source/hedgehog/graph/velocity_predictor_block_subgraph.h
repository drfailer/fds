#ifndef VELOCITY_PREDICTOR_BLOCK_SUBGRAPH_H
#define VELOCITY_PREDICTOR_BLOCK_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/mesh_block_data.h"
#include "../state/mesh_block_state.h"
#include "../fds_fortran_interface.h"

/// Block kernel task for velocity predictor.
/// Processes a K-range sub-block of a single mesh.
class VelocityPredictorBlockKernelTask
    : public hh::AbstractTask<1, MeshBlockData, MeshBlockData> {
public:
    explicit VelocityPredictorBlockKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshBlockData, MeshBlockData>(
              "VelPredBlockKernel", numThreads) {}

    void execute(std::shared_ptr<MeshBlockData> block) override {
        fds_velocity_predictor_block_kernel(block->nm, block->dt,
                                            block->k1, block->k2);
        this->addResult(block);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshBlockData, MeshBlockData>>
    copy() override {
        return std::make_shared<VelocityPredictorBlockKernelTask>(
            this->numberThreads());
    }
};

/// Merged post-reassembly task for velocity predictor.
/// Combines CC_PROJECT_VELOCITY fix + WALL_VELOCITY_NO_GRADH fix + CHECK_STABILITY
/// into a single task to eliminate 2 intermediate queues.
/// CC_PROJECT_VELOCITY and WALL_VELOCITY_NO_GRADH are no-ops for non-CC_IBM / FFT cases.
class VelPredPostReassembleTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit VelPredPostReassembleTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "VelPredPostReassemble", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_cc_project_velocity_kernel(data->nm, data->dt,
                                        /*store=*/0, /*predictor=*/1);
        fds_wall_velocity_no_gradh_kernel(data->nm, data->dt,
                                           /*store=*/0, /*predictor=*/1);
        fds_check_stability_kernel_only(data->nm, data->t + data->dt, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<VelPredPostReassembleTask>(this->numberThreads());
    }
};

/// Build the velocity predictor sub-graph with block decomposition.
///
/// Pipeline:
///   MeshData -> Decompose -> VelPredBlockKernel(parallel) -> Reassemble
///            -> VelPredPostReassemble(CCProjectVel + WallVel + CheckStability)
///            -> MeshData
///
/// @param kernelThreads Number of threads for parallel tasks
/// @param numBlocks Target number of blocks per mesh
inline auto buildVelocityPredictorBlockSubgraph(size_t kernelThreads, int numBlocks) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData, MeshData>>("VelocityPredictorBlock");

    auto decomposeSM = std::make_shared<hh::StateManager<1, MeshData, MeshBlockData>>(
        std::make_shared<MeshBlockDecomposeState>(numBlocks), "VelPredDecompose");
    auto blockKernel = std::make_shared<VelocityPredictorBlockKernelTask>(kernelThreads);
    auto reassembleSM = std::make_shared<hh::StateManager<1, MeshBlockData, MeshData>>(
        std::make_shared<MeshBlockReassembleState>(), "VelPredReassemble");
    auto postReassemble = std::make_shared<VelPredPostReassembleTask>(kernelThreads);

    subgraph->inputs(decomposeSM);
    subgraph->edges(decomposeSM, blockKernel);
    subgraph->edges(blockKernel, reassembleSM);
    subgraph->edges(reassembleSM, postReassemble);
    subgraph->outputs(postReassemble);

    return subgraph;
}

#endif // VELOCITY_PREDICTOR_BLOCK_SUBGRAPH_H
