#ifndef EXCHANGE_FANOUT_TASK_H
#define EXCHANGE_FANOUT_TASK_H

#include <hedgehog/hedgehog.h>
#include <hedgehog_comm.h>
#include "../data/mesh_data.h"
#include "../data/exchange_flux_mesh_data.h"
#include "../tool/mesh_dependency_graph.h"

/// Parallel fan-out task that pushes slab data into ExchangeFluxMeshData objects.
///
/// For each send target of mesh NM, allocates an ExchangeFluxMeshData from the
/// shared pool, copies slab data from Fortran arrays, and emits it downstream
/// to the CommunicatorTask (or WriteBufferTask if no comm service).
///
/// MeshData pending tokens flow to the gate via graph-level input broadcast
/// (the graph input connects to both this task and the gate).
///
/// @tparam Strategy Exchange strategy trait (e.g. FluxExchangeStrategy)
template <typename Strategy>
class ExchangeFanOutTask
    : public hh::AbstractTask<1, MeshData, ExchangeFluxMeshData<Strategy>> {
    using EFD = ExchangeFluxMeshData<Strategy>;
    using Pool = hh::comm::tool::MemoryPool<EFD>;
public:
    ExchangeFanOutTask(size_t numThreads,
                       std::shared_ptr<MeshDependencyGraph> depGraph,
                       std::shared_ptr<Pool> pool)
        : hh::AbstractTask<1, MeshData, EFD>(
              "ExchangeFanOut", numThreads),
          depGraph_(std::move(depGraph)),
          pool_(std::move(pool)) {}

    void execute(std::shared_ptr<MeshData> md) override {
        int nm = md->nm;
        for (int nom : depGraph_->sendTargets(nm)) {
            auto efd = pool_->template allocate<EFD>(
                hh::comm::tool::MemoryManagerAllocateMode::Wait);
            efd->reset(nm, nom);
            efd->push();
            this->addResult(efd);
        }
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, EFD>> copy() override {
        return std::make_shared<ExchangeFanOutTask>(
            this->numberThreads(), depGraph_, pool_);
    }

private:
    std::shared_ptr<MeshDependencyGraph> depGraph_;
    std::shared_ptr<Pool> pool_;
};

#endif // EXCHANGE_FANOUT_TASK_H
