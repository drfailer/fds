#ifndef VELOCITY_CORRECTOR_BLOCK_SUBGRAPH_H
#define VELOCITY_CORRECTOR_BLOCK_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/mesh_block_data.h"
#include "../state/mesh_block_state.h"
#include "../fds_fortran_interface.h"

/// Mesh-level CC_PROJECT_VELOCITY store for CC_IBM (runs before block decomposition).
/// Saves current projected velocities so the corrector fix can compute the average.
/// No-op for non-CC_IBM (checked in Fortran C wrapper).
class CCProjectVelocityCorrStoreTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit CCProjectVelocityCorrStoreTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "CCProjectVel_CorrStore", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_cc_project_velocity_kernel(data->nm, data->dt,
                                        /*store=*/1, /*predictor=*/0);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<CCProjectVelocityCorrStoreTask>(this->numberThreads());
    }
};

/// Mesh-level CC_PROJECT_VELOCITY fix for CC_IBM (runs after block reassembly).
/// Applies projected velocity correction after corrector kernel.
/// No-op for non-CC_IBM (checked in Fortran C wrapper).
class CCProjectVelocityCorrFixTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit CCProjectVelocityCorrFixTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "CCProjectVel_CorrFix", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_cc_project_velocity_kernel(data->nm, data->dt,
                                        /*store=*/0, /*predictor=*/0);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<CCProjectVelocityCorrFixTask>(this->numberThreads());
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

/// Mesh-level wall velocity store for sparse solvers (runs before block decomposition).
/// Saves current wall velocities so the corrector fix can compute the average.
/// No-op for FFT solver (handled inside the Fortran wrapper).
class WallVelNoGradHCorrStoreTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit WallVelNoGradHCorrStoreTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "WallVelNoGradH_CorrStore", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_wall_velocity_no_gradh_kernel(data->nm, data->dt,
                                           /*store=*/1, /*predictor=*/0);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<WallVelNoGradHCorrStoreTask>(this->numberThreads());
    }
};

/// Mesh-level wall velocity fix for sparse solvers (runs after block reassembly).
/// Corrects wall-adjacent velocities using stored values + force terms (no grad H).
/// No-op for FFT solver (handled inside the Fortran wrapper).
class WallVelNoGradHCorrFixTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    explicit WallVelNoGradHCorrFixTask(size_t numThreads)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "WallVelNoGradH_CorrFix", numThreads) {}

    void execute(std::shared_ptr<MeshData> data) override {
        fds_wall_velocity_no_gradh_kernel(data->nm, data->dt,
                                           /*store=*/0, /*predictor=*/0);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>>
    copy() override {
        return std::make_shared<WallVelNoGradHCorrFixTask>(this->numberThreads());
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
///   MeshData -> CCProjectVelStore (no-op for non-CC_IBM)
///            -> WallVelStore (no-op for FFT, stores wall vels for sparse solvers)
///            -> Decompose -> VelCorrBlockKernel(parallel) -> Reassemble
///            -> CCProjectVelFix (no-op for non-CC_IBM)
///            -> WallVelFix (no-op for FFT, fixes wall vels for sparse solvers)
///            -> CheckDivKernel(parallel) -> MeshData
///
/// The velocity corrector kernel (pure I,J,K loops) is block-decomposed along K.
/// The check divergence kernel (reduction to mesh scalars) runs at mesh level after reassembly.
///
/// @param kernelThreads Number of threads for parallel tasks
/// @param numBlocks Target number of blocks per mesh (default: kernelThreads)
inline auto buildVelocityCorrectorBlockSubgraph(size_t kernelThreads, int numBlocks) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData, MeshData>>("VelocityCorrectorBlock");

    auto ccProjectVelStore = std::make_shared<CCProjectVelocityCorrStoreTask>(kernelThreads);
    auto wallVelStore = std::make_shared<WallVelNoGradHCorrStoreTask>(kernelThreads);
    auto decomposeSM = std::make_shared<hh::StateManager<1, MeshData, MeshBlockData>>(
        std::make_shared<MeshBlockDecomposeState>(numBlocks), "VelCorrDecompose");
    auto blockKernel = std::make_shared<VelocityCorrectorBlockKernelTask>(kernelThreads);
    auto reassembleSM = std::make_shared<hh::StateManager<1, MeshBlockData, MeshData>>(
        std::make_shared<MeshBlockReassembleState>(), "VelCorrReassemble");
    auto ccProjectVelFix = std::make_shared<CCProjectVelocityCorrFixTask>(kernelThreads);
    auto wallVelFix = std::make_shared<WallVelNoGradHCorrFixTask>(kernelThreads);
    auto checkDivTask = std::make_shared<CheckDivergenceKernelTask>(kernelThreads);

    subgraph->inputs(ccProjectVelStore);
    subgraph->edges(ccProjectVelStore, wallVelStore);
    subgraph->edges(wallVelStore, decomposeSM);
    subgraph->edges(decomposeSM, blockKernel);
    subgraph->edges(blockKernel, reassembleSM);
    subgraph->edges(reassembleSM, ccProjectVelFix);
    subgraph->edges(ccProjectVelFix, wallVelFix);
    subgraph->edges(wallVelFix, checkDivTask);
    subgraph->outputs(checkDivTask);

    return subgraph;
}

#endif // VELOCITY_CORRECTOR_BLOCK_SUBGRAPH_H
