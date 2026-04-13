#ifndef EXCHANGE_WRITE_BUFFER_TASK_H
#define EXCHANGE_WRITE_BUFFER_TASK_H

#include <hedgehog/hedgehog.h>
#include <hedgehog_comm.h>
#include <cstring>
#include <array>
#include <unordered_map>
#include "../data/exchange_mesh_data.h"
#include "../tool/exchange_buffer.h"

/// Task that receives ExchangeMeshData from the CommunicatorTask (or
/// directly from FanOutTask), writes the slab data into the shared
/// ExchangeBuffer for the appropriate exchange code and round slot,
/// releases the ExchangeMeshData back to the pool, and emits an
/// ExchangeDepSignal to the gate.
///
/// Double-buffered: uses roundId % 2 to select the buffer slot, so
/// concurrent exchange rounds never share storage.
class ExchangeWriteBufferTask
    : public hh::AbstractTask<1, ExchangeMeshData, ExchangeDepSignal> {
    using Pool = hh::comm::tool::MemoryPool<ExchangeMeshData>;
    using BufferMap = std::unordered_map<int, std::array<std::shared_ptr<ExchangeBuffer>, 2>>;
public:
    ExchangeWriteBufferTask(
        size_t numThreads,
        std::shared_ptr<BufferMap> buffers,
        std::shared_ptr<Pool> pool)
        : hh::AbstractTask<1, ExchangeMeshData, ExchangeDepSignal>(
              "ExchangeWriteBuffer", numThreads),
          buffers_(std::move(buffers)),
          pool_(std::move(pool)) {}

    void execute(std::shared_ptr<ExchangeMeshData> emd) override {
        int nm    = emd->header_.sourceNm;
        int nom   = emd->header_.destNom;
        int sz    = emd->header_.slabSize;
        int code  = emd->header_.exchangeCode;
        int round = emd->header_.roundId;
        if (sz > 0) {
            auto &buf = buffers_->at(code)[round % 2];
            std::memcpy(buf->buffer(nm, nom), emd->buffer_.data(),
                        static_cast<size_t>(sz) * sizeof(double));
        }
        pool_->template release<ExchangeMeshData>(std::move(emd));
        this->addResult(std::make_shared<ExchangeDepSignal>(nm, nom, round));
    }

    std::shared_ptr<hh::AbstractTask<1, ExchangeMeshData, ExchangeDepSignal>> copy() override {
        return std::make_shared<ExchangeWriteBufferTask>(
            this->numberThreads(), buffers_, pool_);
    }

private:
    std::shared_ptr<BufferMap> buffers_;
    std::shared_ptr<Pool> pool_;
};

#endif // EXCHANGE_WRITE_BUFFER_TASK_H
