#ifndef EXCHANGE_FANOUT_TASK_H
#define EXCHANGE_FANOUT_TASK_H

#include <hedgehog/hedgehog.h>
#include <hedgehog_comm.h>
#include "../data/mesh_data.h"
#include "../data/exchange_mesh_data.h"
#include "../tool/mesh_dependency_graph.h"

/// Parallel fan-out task that pushes slab data into ExchangeMeshData objects.
///
/// For each send target of mesh NM, allocates an ExchangeMeshData from the
/// shared pool, copies slab data from Fortran arrays, and emits it downstream
/// to the CommunicatorTask (or WriteBufferTask if no comm service).
///
/// The exchange code is read from md->exchangeCode at runtime.
///
/// MeshData pending tokens flow to the gate via graph-level input broadcast
/// (the graph input connects to both this task and the gate).
template<MeshState S = MeshState::Default>
class ExchangeFanOutTask
    : public hh::AbstractTask<1, MeshData<S>, ExchangeMeshData> {
    using Pool = hh::comm::tool::MemoryPool<ExchangeMeshData>;
public:
    ExchangeFanOutTask(size_t numThreads,
                       std::shared_ptr<MeshDependencyGraph> depGraph,
                       std::shared_ptr<Pool> pool)
        : hh::AbstractTask<1, MeshData<S>, ExchangeMeshData>(
              "ExchangeFanOut", numThreads),
          depGraph_(std::move(depGraph)),
          pool_(std::move(pool)) {}

    void execute(std::shared_ptr<MeshData<S>> md) override {
        int nm = md->nm;
        int code = md->exchangeCode;
        int round = md->exchangeRound;
        for (int nom : depGraph_->sendTargets(nm)) {
            auto emd = pool_->template allocate<ExchangeMeshData>(
                hh::comm::tool::MemoryManagerAllocateMode::Wait);
            emd->reset(nm, nom, code);
            emd->header_.roundId = round;
            emd->push();
            this->addResult(emd);
        }
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<S>, ExchangeMeshData>> copy() override {
        return std::make_shared<ExchangeFanOutTask<S>>(
            this->numberThreads(), depGraph_, pool_);
    }

private:
    std::shared_ptr<MeshDependencyGraph> depGraph_;
    std::shared_ptr<Pool> pool_;
};

#endif // EXCHANGE_FANOUT_TASK_H
