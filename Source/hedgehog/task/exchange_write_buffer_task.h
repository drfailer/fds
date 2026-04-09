#ifndef EXCHANGE_WRITE_BUFFER_TASK_H
#define EXCHANGE_WRITE_BUFFER_TASK_H

#include <hedgehog/hedgehog.h>
#include <hedgehog_comm.h>
#include <cstring>
#include "../data/exchange_flux_mesh_data.h"
#include "../tool/exchange_buffer.h"

/// Task that receives ExchangeFluxMeshData from the CommunicatorTask (or
/// directly from FanOutTask), writes the slab data into the shared
/// ExchangeBuffer, releases the ExchangeFluxMeshData back to the pool,
/// and emits an ExchangeDepSignal to the gate.
///
/// @tparam Strategy Exchange strategy trait (e.g. FluxExchangeStrategy)
template <typename Strategy>
class ExchangeWriteBufferTask
    : public hh::AbstractTask<1, ExchangeFluxMeshData<Strategy>, ExchangeDepSignal> {
    using EFD = ExchangeFluxMeshData<Strategy>;
    using Pool = hh::comm::tool::MemoryPool<EFD>;
public:
    ExchangeWriteBufferTask(size_t numThreads,
                            std::shared_ptr<ExchangeBuffer<Strategy>> buffer,
                            std::shared_ptr<Pool> pool)
        : hh::AbstractTask<1, EFD, ExchangeDepSignal>(
              "ExchangeWriteBuffer", numThreads),
          buffer_(std::move(buffer)),
          pool_(std::move(pool)) {}

    void execute(std::shared_ptr<EFD> efd) override {
        int nm  = efd->header_.sourceNm;
        int nom = efd->header_.destNom;
        int sz  = efd->header_.slabSize;
        if (sz > 0) {
            std::memcpy(buffer_->buffer(nm, nom), efd->buffer_.data(),
                        static_cast<size_t>(sz) * sizeof(double));
        }
        pool_->template release<EFD>(std::move(efd));
        this->addResult(std::make_shared<ExchangeDepSignal>(nm, nom));
    }

    std::shared_ptr<hh::AbstractTask<1, EFD, ExchangeDepSignal>> copy() override {
        return std::make_shared<ExchangeWriteBufferTask>(
            this->numberThreads(), buffer_, pool_);
    }

private:
    std::shared_ptr<ExchangeBuffer<Strategy>> buffer_;
    std::shared_ptr<Pool> pool_;
};

#endif // EXCHANGE_WRITE_BUFFER_TASK_H
