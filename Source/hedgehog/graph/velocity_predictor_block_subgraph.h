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

/// Mesh-level CFL/VN stability check task (runs after block reassembly).
class CheckStabilityKernelTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit CheckStabilityKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "CheckStabilityKernel", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_check_stability_kernel_only(data->nm, data->t + data->dt, data->dt);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<CheckStabilityKernelTask>(
            this->numberThreads());
    }
};

/// Build the velocity predictor sub-graph with block decomposition.
///
/// Pipeline:
///   MeshData -> Decompose -> VelPredBlockKernel(parallel) -> Reassemble
///            -> CheckStabilityKernel(parallel, if !skipCFL) -> MeshData
///
/// Note: This sub-graph does NOT call WALL_VELOCITY_NO_GRADH, so it must
/// not be used for sparse pressure solvers (ULMAT/GLMAT/UGLMAT) which
/// require that call between the kernel and CHECK_STABILITY.
/// For those solvers, use VelocityPredictorFullTask instead.
///
/// @param kernelThreads Number of threads for parallel tasks
/// @param numBlocks Target number of blocks per mesh
/// @param skipCFL If true, skip CHECK_STABILITY_KERNEL (CC_IBM path)
inline auto buildVelocityPredictorBlockSubgraph(size_t kernelThreads, int numBlocks,
                                                 bool skipCFL) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData, MeshData>>("VelocityPredictorBlock");

    auto decomposeSM = std::make_shared<hh::StateManager<1, MeshData, MeshBlockData>>(
        std::make_shared<MeshBlockDecomposeState>(numBlocks), "VelPredDecompose");
    auto blockKernel = std::make_shared<VelocityPredictorBlockKernelTask>(kernelThreads);
    auto reassembleSM = std::make_shared<hh::StateManager<1, MeshBlockData, MeshData>>(
        std::make_shared<MeshBlockReassembleState>(), "VelPredReassemble");

    subgraph->inputs(decomposeSM);
    subgraph->edges(decomposeSM, blockKernel);
    subgraph->edges(blockKernel, reassembleSM);

    if (skipCFL) {
        subgraph->outputs(reassembleSM);
    } else {
        auto checkStabilityTask = std::make_shared<CheckStabilityKernelTask>(kernelThreads);
        subgraph->edges(reassembleSM, checkStabilityTask);
        subgraph->outputs(checkStabilityTask);
    }

    return subgraph;
}

#endif // VELOCITY_PREDICTOR_BLOCK_SUBGRAPH_H
