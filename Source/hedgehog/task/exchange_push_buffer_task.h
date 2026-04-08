#ifndef EXCHANGE_PUSH_BUFFER_TASK_H
#define EXCHANGE_PUSH_BUFFER_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../tool/mesh_dependency_graph.h"
#include "../tool/exchange_buffer.h"

/// Parallel task that pushes source mesh arrays into the exchange buffer.
///
/// For each send target of mesh NM, copies the exchange-relevant arrays
/// (determined by Strategy) from MESHES(NM) into the shared ExchangeBuffer.
///
/// Thread safety: multiple threads push different source meshes (NM)
/// concurrently.  Each writes to non-overlapping buffer entries.
///
/// @tparam Strategy Exchange strategy trait (e.g. FluxExchangeStrategy)
template <typename Strategy>
class ExchangePushBufferTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    ExchangePushBufferTask(size_t numThreads,
                           std::shared_ptr<MeshDependencyGraph> depGraph,
                           std::shared_ptr<ExchangeBuffer<Strategy>> buffer)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "ExchangePush", numThreads),
          depGraph_(std::move(depGraph)),
          buffer_(std::move(buffer)) {}

    void execute(std::shared_ptr<MeshData> md) override {
        int nm = md->nm;
        for (int nom : depGraph_->sendTargets(nm)) {
            int sz = buffer_->bufferSize(nm, nom);
            if (sz > 0) {
                Strategy::pushToBuffer(nm, nom, buffer_->buffer(nm, nom), sz);
            }
        }
        this->addResult(md);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>> copy() override {
        return std::make_shared<ExchangePushBufferTask>(
            this->numberThreads(), depGraph_, buffer_);
    }

private:
    std::shared_ptr<MeshDependencyGraph> depGraph_;
    std::shared_ptr<ExchangeBuffer<Strategy>> buffer_;
};

#endif // EXCHANGE_PUSH_BUFFER_TASK_H
