#ifndef VISC_DENSITY_BLOCK_SUBGRAPH_H
#define VISC_DENSITY_BLOCK_SUBGRAPH_H

#include <hedgehog/hedgehog.h>
#include <algorithm>
#include <memory>
#include <unordered_map>
#include "../data/mesh_data.h"
#include "../data/mesh_block_data.h"
#include "../fds_fortran_interface.h"
#include "compute_viscosity_block_subgraph.h"
#include "density_block_subgraph.h"

/// Mid-state bridging viscosity and density block kernels.
///
/// Reassembles visc K-blocks per-mesh, then runs three sequential operations:
///   1. ViscPostBlock (wall loops + corner mirroring)
///   2. MassFD (diffusion terms)
///   3. DensityPreprocessing (work arrays, settling velocity, wall loop)
/// Then immediately re-decomposes for the density kernel.
///
/// Dispatches per-mesh (does NOT wait for all N meshes), so density blocks for
/// mesh 1 can execute in parallel with visc blocks for mesh 2.
class ViscDensityMidState
    : public hh::AbstractState<1, MeshBlockData, MeshBlockData> {
public:
    explicit ViscDensityMidState(int numBlocks)
        : numBlocks_(std::max(1, numBlocks)) {}

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

            // Sequential per-mesh operations
            fds_compute_viscosity_post_block(md->nm, md->phase);
            fds_mass_finite_differences_kernel(md->nm);
            fds_density_block_preprocessing(md->nm, md->t, md->dt);

            // Re-decompose for density kernel (same K-block scheme)
            int kbar = fds_get_kbar(md->nm);
            int bs = std::max(1, (kbar + numBlocks_ - 1) / numBlocks_);
            int total = (kbar + bs - 1) / bs;
            for (int b = 0; b < total; ++b) {
                int k1 = b * bs + 1;
                int k2 = std::min((b + 1) * bs, kbar);
                this->addResult(std::make_shared<MeshBlockData>(
                    md->nm, k1, k2, md->t, md->dt, md->phase, total, md));
            }

            entries_.erase(nm);
        }
    }

private:
    int numBlocks_;
    struct Entry {
        int count = 0;
        int expected = 0;
        std::shared_ptr<MeshData> meshData;
    };
    std::unordered_map<int, Entry> entries_;
};

/// Build the merged viscosity + density sub-graph with block decomposition.
///
/// Pipeline:
///   MeshData -> ViscOrch(CC cutface, K-decompose) ->
///   ViscBlockKernel(parallel: MU_DNS, strain rate, turb MU, KRES) ->
///   MidState(reassemble, ViscPostBlock + MassFD + DensityPreproc, re-decompose) ->
///   DensityBlockKernel(parallel: species density, M_DOT_PPP, RHOS/RHO) ->
///   DensityCollector(STORE_FLUX, CHECK_MASS, ZZ/=RHO, CLIP, PBAR, RSUM, TMP)
///   -> MeshData
///
/// Eliminates 2 nodes and 2 queues vs separate ViscBlock + MassFD + DensityBlock.
inline auto buildViscDensityBlockSubgraph(int nmeshes, size_t kernelThreads,
                                           int numBlocks) {
    auto subgraph = std::make_shared<hh::Graph<1, MeshData, MeshData>>("ViscDensityBlock");

    // Reuse existing orchestrator, kernel tasks, and collector classes
    auto orchestratorSM = std::make_shared<hh::StateManager<1, MeshData, MeshBlockData>>(
        std::make_shared<ComputeViscosityBlockOrchestrator>(nmeshes, numBlocks),
        "ViscDensityOrch");
    auto viscKernel = std::make_shared<ComputeViscosityBlockKernelTask>(kernelThreads);
    auto midStateSM = std::make_shared<hh::StateManager<1, MeshBlockData, MeshBlockData>>(
        std::make_shared<ViscDensityMidState>(numBlocks), "ViscDensityMid");
    auto densityKernel = std::make_shared<DensityBlockKernelTask>(kernelThreads);
    auto collectorSM = std::make_shared<hh::StateManager<1, MeshBlockData, MeshData>>(
        std::make_shared<DensityBlockCollector>(), "ViscDensityCollector");

    subgraph->inputs(orchestratorSM);
    subgraph->edges(orchestratorSM, viscKernel);
    subgraph->edges(viscKernel, midStateSM);
    subgraph->edges(midStateSM, densityKernel);
    subgraph->edges(densityKernel, collectorSM);
    subgraph->outputs(collectorSM);

    return subgraph;
}

#endif // VISC_DENSITY_BLOCK_SUBGRAPH_H
