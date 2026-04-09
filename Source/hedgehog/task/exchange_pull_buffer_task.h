#ifndef EXCHANGE_PULL_BUFFER_TASK_H
#define EXCHANGE_PULL_BUFFER_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../tool/exchange_buffer.h"

/// Parallel task that pulls exchange data from the buffer into OMESH arrays.
///
/// For each source mesh NM that sends to this mesh NOM (same-rank and
/// cross-rank), copies the buffered data into MESHES(NOM)%OMESH(NM)
/// using the Strategy's pullFromBuffer.
///
/// Thread safety: multiple threads pull for different destination meshes (NOM)
/// concurrently.  Each writes to non-overlapping OMESH entries.
///
/// @tparam Strategy Exchange strategy trait (e.g. FluxExchangeStrategy)
template <typename Strategy>
class ExchangePullBufferTask
    : public hh::AbstractTask<1, MeshData, MeshData> {
public:
    ExchangePullBufferTask(size_t numThreads,
                           std::shared_ptr<ExchangeBuffer<Strategy>> buffer)
        : hh::AbstractTask<1, MeshData, MeshData>(
              "ExchangePull", numThreads),
          buffer_(std::move(buffer)) {}

    void execute(std::shared_ptr<MeshData> md) override {
        int nom = md->nm;  // this mesh is the receiver
        for (int nm : buffer_->recvSources(nom)) {
            int sz = buffer_->bufferSize(nm, nom);
            if (sz > 0) {
                Strategy::pullFromBufferRecv(nom, nm, buffer_->buffer(nm, nom), sz);
            }
        }
        this->addResult(md);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData, MeshData>> copy() override {
        return std::make_shared<ExchangePullBufferTask>(
            this->numberThreads(), buffer_);
    }

private:
    std::shared_ptr<ExchangeBuffer<Strategy>> buffer_;
};

#endif // EXCHANGE_PULL_BUFFER_TASK_H
