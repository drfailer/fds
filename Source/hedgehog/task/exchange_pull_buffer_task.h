#ifndef EXCHANGE_PULL_BUFFER_TASK_H
#define EXCHANGE_PULL_BUFFER_TASK_H

#include <hedgehog/hedgehog.h>
#include <array>
#include <unordered_map>
#include "../data/mesh_data.h"
#include "../tool/exchange_buffer.h"
#include "../fds_fortran_interface.h"

/// Parallel task that pulls exchange data from the buffer into OMESH arrays.
///
/// Double-buffered: uses exchangeRound % 2 to select the buffer slot,
/// matching the slot used by WriteBufferTask for the same round.
///
/// Thread safety: multiple threads pull for different destination meshes (NOM)
/// concurrently.  Each writes to non-overlapping OMESH entries.
template<MeshState S = MeshState::Default>
class ExchangePullBufferTask
    : public hh::AbstractTask<1, MeshData<S>, MeshData<S>> {
    using BufferMap = std::unordered_map<int, std::array<std::shared_ptr<ExchangeBuffer>, 2>>;
public:
    ExchangePullBufferTask(
        size_t numThreads,
        std::shared_ptr<BufferMap> buffers)
        : hh::AbstractTask<1, MeshData<S>, MeshData<S>>(
              "ExchangePull", numThreads),
          buffers_(std::move(buffers)) {}

    void execute(std::shared_ptr<MeshData<S>> md) override {
        int nom = md->nm;  // this mesh is the receiver
        int code = md->exchangeCode;
        int slot = md->exchangeRound % 2;
        auto &buf = buffers_->at(code)[slot];
        for (int nm : buf->recvSources(nom)) {
            int sz = buf->bufferSize(nm, nom);
            if (sz > 0) {
                fds_exchange_pull_slab_recv(code, nom, nm, buf->buffer(nm, nom), sz);
            }
        }
        this->addResult(md);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<S>, MeshData<S>>> copy() override {
        return std::make_shared<ExchangePullBufferTask<S>>(
            this->numberThreads(), buffers_);
    }

private:
    std::shared_ptr<BufferMap> buffers_;
};

#endif // EXCHANGE_PULL_BUFFER_TASK_H
