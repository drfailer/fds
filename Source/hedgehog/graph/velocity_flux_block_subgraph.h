#ifndef VELOCITY_FLUX_BLOCK_SUBGRAPH_H
#define VELOCITY_FLUX_BLOCK_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <memory>
#include "../data/mesh_data.h"
#include "../data/mesh_block_data.h"
#include "../state/mesh_block_state.h"
#include "../fds_fortran_interface.h"

/// Block kernel task for velocity flux (DivSetup).
/// Processes a K-range sub-block of a single mesh: vorticity, FVX, FVY, FVZ, DIRECT_FORCE.
class VelocityFluxBlockKernelTask
    : public hh::AbstractTask<1, MeshBlockData, MeshBlockData> {
public:
    explicit VelocityFluxBlockKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshBlockData, MeshBlockData>(
              "VelFluxBlockKernel", numThreads) {}

    void execute(std::shared_ptr<MeshBlockData> block) override {
        fds_velocity_flux_block_kernel(block->nm, block->t, block->dt,
                                       block->phase, block->k1, block->k2);
        this->addResult(block);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshBlockData, MeshBlockData>>
    copy() override {
        return std::make_shared<VelocityFluxBlockKernelTask>(
            this->numberThreads());
    }
};

/// Orchestrator state for velocity flux block decomposition.
/// Collects N MeshData tokens, runs sequential pre-processing per mesh,
/// then decomposes each mesh into K-blocks for parallel execution.
class VelocityFluxBlockOrchestrator
    : public hh::AbstractState<1, MeshData, MeshBlockData> {
public:
    VelocityFluxBlockOrchestrator(int nmeshes, int numBlocks)
        : hh::AbstractState<1, MeshData, MeshBlockData>(),
          nmeshes_(nmeshes), numBlocks_(std::max(1, numBlocks)) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // Sequential pre-processing per mesh
            for (auto &md : collected_) {
                // CC_IBM: CC_VELOCITY_BC must run before velocity flux
                if (fds_is_cc_ibm())
                    fds_cc_velocity_bc(md->t, md->nm, md->phase);
                fds_set_baroclinic_false(md->nm);
                fds_viscosity_bc_kernel(md->nm, md->phase);
                // CC_IBM: set cutface velocities before block kernels
                if (fds_is_cc_ibm())
                    fds_cutface_velocities(md->nm, md->phase, 1);
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

/// Collector state for velocity flux block decomposition.
/// Reassembles blocks back into MeshData, runs sequential post-processing.
class VelocityFluxBlockCollector
    : public hh::AbstractState<1, MeshBlockData, MeshData> {
public:
    VelocityFluxBlockCollector() = default;

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
            // CC_IBM post-processing: reset cutface velocities + CC gravity corrections
            if (fds_is_cc_ibm())
                fds_cc_velocity_flux_post(md->nm, md->t, md->dt, md->phase);
            // Post-processing: agglomeration (corrector only)
            if (md->phase)
                fds_agglomeration(md->dt, md->nm);
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

/// Build the velocity flux (DivSetup) sub-graph with block decomposition.
///
/// Pipeline:
///   MeshData -> Orchestrator(baroclinic+viscBC, decompose) ->
///   VelFluxBlockKernel(parallel) -> Collector(agglomeration) -> MeshData
///
/// @param nmeshes Number of meshes
/// @param kernelThreads Number of threads for parallel tasks
/// @param numBlocks Target number of blocks per mesh
inline auto buildVelocityFluxBlockSubgraph(int nmeshes, size_t kernelThreads,
                                            int numBlocks) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData, MeshData>>("VelocityFluxBlock");

    auto orchestratorSM = std::make_shared<hh::StateManager<1, MeshData, MeshBlockData>>(
        std::make_shared<VelocityFluxBlockOrchestrator>(nmeshes, numBlocks),
        "VelFluxOrch");
    auto blockKernel = std::make_shared<VelocityFluxBlockKernelTask>(kernelThreads);
    auto collectorSM = std::make_shared<hh::StateManager<1, MeshBlockData, MeshData>>(
        std::make_shared<VelocityFluxBlockCollector>(), "VelFluxCollector");

    subgraph->inputs(orchestratorSM);
    subgraph->edges(orchestratorSM, blockKernel);
    subgraph->edges(blockKernel, collectorSM);
    subgraph->outputs(collectorSM);

    return subgraph;
}

#endif // VELOCITY_FLUX_BLOCK_SUBGRAPH_H
