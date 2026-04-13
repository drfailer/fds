// UNUSED — Block decomposition disabled. Kept for reference.
#ifndef VELOCITY_CORRECTOR_BLOCK_SUBGRAPH_H
#define VELOCITY_CORRECTOR_BLOCK_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/mesh_block_data.h"
#include "../state/mesh_block_state.h"
#include "../fds_fortran_interface.h"

/// Merged pre-decompose kernel for velocity corrector.
/// CC_PROJECT_VELOCITY store + WALL_VELOCITY_NO_GRADH store.
class VelCorrPreKernelTask
    : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit VelCorrPreKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>(
              "VelCorrPreKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_cc_project_velocity_kernel(data->nm, data->dt, 1, 0);
        fds_wall_velocity_no_gradh_kernel(data->nm, data->dt, 1, 0);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>>
    copy() override {
        return std::make_shared<VelCorrPreKernelTask>(this->numberThreads());
    }
};

/// Block kernel task for velocity corrector.
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

/// Merged post-reassembly kernel for velocity corrector.
/// CC_PROJECT_VELOCITY fix + WALL_VELOCITY_NO_GRADH fix + CHECK_DIVERGENCE.
class VelCorrPostKernelTask
    : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit VelCorrPostKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>(
              "VelCorrPostKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_cc_project_velocity_kernel(data->nm, data->dt, 0, 0);
        fds_wall_velocity_no_gradh_kernel(data->nm, data->dt, 0, 0);
        fds_check_divergence_kernel(data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>>
    copy() override {
        return std::make_shared<VelCorrPostKernelTask>(this->numberThreads());
    }
};

/// Build the velocity corrector sub-graph with block decomposition.
///
/// Pipeline:
///   MeshData<> -> VelCorrPreKernel(CCProjectVelStore + WallVelStore)
///            -> Decompose -> VelCorrBlockKernel(parallel) -> Reassemble
///            -> VelCorrPostKernel(CCProjectVelFix + WallVelFix + CheckDiv)
///            -> MeshData<>
inline auto buildVelocityCorrectorBlockSubgraph(size_t blockThreads,
                                                  int numBlocks, int nmeshes) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData<>, MeshData<>>>("VelocityCorrectorBlock");

    auto preKernel = std::make_shared<VelCorrPreKernelTask>(static_cast<size_t>(nmeshes));
    auto decomposeTask = std::make_shared<MeshBlockDecomposeTask>(numBlocks, "VelCorrDecompose");
    auto blockKernel = std::make_shared<VelocityCorrectorBlockKernelTask>(blockThreads);
    auto reassembleTask = std::make_shared<MeshBlockReassembleTask>("VelCorrReassemble");
    auto postKernel = std::make_shared<VelCorrPostKernelTask>(static_cast<size_t>(nmeshes));

    subgraph->inputs(preKernel);
    subgraph->edges(preKernel, decomposeTask);
    subgraph->edges(decomposeTask, blockKernel);
    subgraph->edges(blockKernel, reassembleTask);
    subgraph->edges(reassembleTask, postKernel);
    subgraph->outputs(postKernel);

    return subgraph;
}

#endif // VELOCITY_CORRECTOR_BLOCK_SUBGRAPH_H
