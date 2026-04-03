// UNUSED — Block decomposition disabled. Kept for reference.
#ifndef VELOCITY_PREDICTOR_BLOCK_SUBGRAPH_H
#define VELOCITY_PREDICTOR_BLOCK_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/mesh_block_data.h"
#include "../state/mesh_block_state.h"
#include "../fds_fortran_interface.h"

/// Block kernel task for velocity predictor.
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

/// Merged post-reassembly kernel for velocity predictor.
/// CC_PROJECT_VELOCITY fix + WALL_VELOCITY_NO_GRADH fix + CHECK_STABILITY.
class VelPredPostKernelTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit VelPredPostKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "VelPredPostKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_cc_project_velocity_kernel(data->nm, data->dt, 0, 1);
        fds_wall_velocity_no_gradh_kernel(data->nm, data->dt, 0, 1);
        fds_check_stability_kernel_only(data->nm, data->t + data->dt, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<VelPredPostKernelTask>(this->numberThreads());
    }
};

/// Build the velocity predictor sub-graph with block decomposition.
///
/// Pipeline:
///   MeshData -> Decompose -> VelPredBlockKernel(parallel) -> Reassemble
///            -> VelPredPostKernel(CCProjectVel + WallVel + CheckStability)
///            -> MeshData
inline auto buildVelocityPredictorBlockSubgraph(size_t blockThreads,
                                                  int numBlocks, int nmeshes) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData, MeshData>>("VelocityPredictorBlock");

    auto decomposeTask = std::make_shared<MeshBlockDecomposeTask>(numBlocks, "VelPredDecompose");
    auto blockKernel = std::make_shared<VelocityPredictorBlockKernelTask>(blockThreads);
    auto reassembleTask = std::make_shared<MeshBlockReassembleTask>("VelPredReassemble");
    auto postKernel = std::make_shared<VelPredPostKernelTask>(static_cast<size_t>(nmeshes));

    subgraph->inputs(decomposeTask);
    subgraph->edges(decomposeTask, blockKernel);
    subgraph->edges(blockKernel, reassembleTask);
    subgraph->edges(reassembleTask, postKernel);
    subgraph->outputs(postKernel);

    return subgraph;
}

#endif // VELOCITY_PREDICTOR_BLOCK_SUBGRAPH_H
