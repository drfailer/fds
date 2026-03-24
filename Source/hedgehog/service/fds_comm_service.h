#ifndef FDS_COMM_SERVICE_H
#define FDS_COMM_SERVICE_H

#include <service/comm_service.hpp>
#include <service/request.hpp>
#include <protocol.hpp>
#include <cassert>
#include <cstddef>
#include <mutex>
#include <vector>
#include <mpi.h>

/// CommService implementation that reuses FDS's already-initialized MPI.
///
/// Unlike MPIService, this does NOT call MPI_Init or MPI_Finalize.
/// FDS initializes MPI in its own Fortran initialization sequence;
/// this service wraps those MPI calls for the communicator task.
class FDSMPIService : public hh::comm::CommService {
public:
    explicit FDSMPIService(bool profilingEnabled = false)
        : hh::comm::CommService(profilingEnabled) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        MPI_Comm_size(MPI_COMM_WORLD, &nbProcesses_);
        // Channel 0 = MPI_COMM_WORLD (matches MPIService convention)
        comms_.emplace_back(MPI_COMM_WORLD);
    }

    ~FDSMPIService() override = default;

private:
    struct MPIRequest {
        MPI_Request request;
        MPI_Status  status;
        MPI_Comm    comm;
        int         flag;
    };

    // Header encoding: matches MPIService's HEADER_FIELDS layout
    static constexpr hh::comm::Header::FieldInfo HEADER_FIELDS[]{
        {.offset = 32, .mask = 0b1111111111111111111111111111111111111111111111111000000000000000}, // source
        {.offset = 32, .mask = 0b1111111111111111111111111111111111111111111111111000000000000000}, // channel
        {.offset = 14, .mask = 0b0000000000000000000000000000000000000000000000000100000000000000}, // signal
        {.offset = 11, .mask = 0b0000000000000000000000000000000000000000000000000011100000000000}, // typeid
        {.offset = 0,  .mask = 0b0000000000000000000000000000000000000000000000000000000000000011}, // buffer id
    };

    static int headerToTag(hh::comm::Header const &header) {
        std::uint64_t tag = 0;
        tag |= header.signal   << HEADER_FIELDS[hh::comm::Header::SIGNAL].offset;
        tag |= header.typeId   << HEADER_FIELDS[hh::comm::Header::TYPE_ID].offset;
        tag |= header.bufferId << HEADER_FIELDS[hh::comm::Header::BUFFER_ID].offset;
        assert((tag & HEADER_FIELDS[0].mask) == 0);
        return (int)tag;
    }

    static hh::comm::Header tagToHeader(int tag) {
        assert(tag >= 0);
        hh::comm::Header header;
        header.signal   = (tag & HEADER_FIELDS[hh::comm::Header::SIGNAL].mask) >> HEADER_FIELDS[hh::comm::Header::SIGNAL].offset;
        header.typeId   = (tag & HEADER_FIELDS[hh::comm::Header::TYPE_ID].mask) >> HEADER_FIELDS[hh::comm::Header::TYPE_ID].offset;
        header.bufferId = (tag & HEADER_FIELDS[hh::comm::Header::BUFFER_ID].mask) >> HEADER_FIELDS[hh::comm::Header::BUFFER_ID].offset;
        return header;
    }

public:
    // --- send ---
    void send(hh::comm::Header const &header, hh::comm::rank_t dest,
              hh::comm::Buffer const &buffer) override {
        std::lock_guard<std::mutex> lk(this->mutex());
        int tag = headerToTag(header);
        MPI_Send(buffer.data(), (int)buffer.size(), MPI_BYTE,
                 (int)dest, tag, comms_[header.channel]);
    }

    hh::comm::Request sendAsync(hh::comm::Header const &header, hh::comm::rank_t dest,
                                hh::comm::Buffer const &buffer) override {
        std::lock_guard<std::mutex> lk(this->mutex());
        auto *r = requestPool_.allocate();
        int tag = headerToTag(header);
        r->comm = comms_[header.channel];
        MPI_Isend(buffer.data(), (int)buffer.size(), MPI_BYTE,
                  (int)dest, tag, r->comm, &r->request);
        return r;
    }

    // --- recv ---
    void recv(hh::comm::Header const &header, hh::comm::Buffer const &buffer) override {
        std::lock_guard<std::mutex> lk(this->mutex());
        MPI_Status status;
        int tag = headerToTag(header);
        MPI_Recv(buffer.data(), (int)buffer.size(), MPI_BYTE,
                 (int)header.source, tag, comms_[header.channel], &status);
    }

    hh::comm::Request recvAsync(hh::comm::Header const &header,
                                hh::comm::Buffer const &buffer) override {
        std::lock_guard<std::mutex> lk(this->mutex());
        auto *r = requestPool_.allocate();
        int tag = headerToTag(header);
        r->comm = comms_[header.channel];
        MPI_Irecv(buffer.data(), (int)buffer.size(), MPI_BYTE,
                  (int)header.source, tag, r->comm, &r->request);
        return r;
    }

    void recv(hh::comm::Request probeRequest, hh::comm::Buffer const &buffer) override {
        std::lock_guard<std::mutex> lk(this->mutex());
        auto *r = static_cast<MPIRequest *>(probeRequest);
        MPI_Recv(buffer.data(), (int)buffer.size(), MPI_BYTE,
                 r->status.MPI_SOURCE, r->status.MPI_TAG, r->comm, &r->status);
        requestPool_.release(r);
    }

    hh::comm::Request recvAsync(hh::comm::Request probeRequest,
                                hh::comm::Buffer const &buffer) override {
        std::lock_guard<std::mutex> lk(this->mutex());
        auto *r = static_cast<MPIRequest *>(probeRequest);
        MPI_Irecv(buffer.data(), (int)buffer.size(), MPI_BYTE,
                  r->status.MPI_SOURCE, r->status.MPI_TAG, r->comm, &r->request);
        return probeRequest;
    }

    // --- probe ---
    hh::comm::Request probe(hh::comm::channel_t channel) override {
        return probe(channel, (hh::comm::rank_t)MPI_ANY_SOURCE);
    }

    hh::comm::Request probeAsync(hh::comm::channel_t channel) override {
        return probeAsync(channel, (hh::comm::rank_t)MPI_ANY_SOURCE);
    }

    hh::comm::Request probe(hh::comm::channel_t channel, hh::comm::rank_t source) override {
        std::lock_guard<std::mutex> lk(this->mutex());
        auto *r = requestPool_.allocate();
        r->comm = comms_[channel];
        MPI_Probe((int)source, MPI_ANY_TAG, r->comm, &r->status);
        return r;
    }

    hh::comm::Request probeAsync(hh::comm::channel_t channel, hh::comm::rank_t source) override {
        std::lock_guard<std::mutex> lk(this->mutex());
        auto *r = requestPool_.allocate();
        r->comm = comms_[channel];
        MPI_Iprobe((int)source, MPI_ANY_TAG, r->comm, &r->flag, &r->status);
        return r;
    }

    // --- request management ---
    bool requestCompleted(hh::comm::Request request) override {
        std::lock_guard<std::mutex> lk(this->mutex());
        auto *r = static_cast<MPIRequest *>(request);
        r->flag = 0;
        MPI_Test(&r->request, &r->flag, &r->status);
        return r->flag != 0;
    }

    void requestRelease(hh::comm::Request request) override {
        requestPool_.release(static_cast<MPIRequest *>(request));
    }

    void requestCancel(hh::comm::Request request) override {
        std::lock_guard<std::mutex> lk(this->mutex());
        auto *r = static_cast<MPIRequest *>(request);
        MPI_Cancel(&r->request);
        requestPool_.release(r);
    }

    size_t bufferSize(hh::comm::Request request) override {
        int count = -1;
        auto *r = static_cast<MPIRequest *>(request);
        MPI_Get_count(&r->status, MPI_BYTE, &count);
        assert(count > 0);
        return (size_t)count;
    }

    hh::comm::Header requestHeader(hh::comm::Request request) override {
        auto *r = static_cast<MPIRequest *>(request);
        auto header = tagToHeader(r->status.MPI_TAG);
        header.channel = (hh::comm::channel_t)r->comm;
        header.source = r->status.MPI_SOURCE;
        return header;
    }

    bool probeSuccess(hh::comm::Request request) override {
        auto *r = static_cast<MPIRequest *>(request);
        return r->flag != 0;
    }

    // --- synchronization ---
    void barrier(hh::comm::channel_t channel = 0) override {
        MPI_Barrier(comms_[channel]);
    }

    // --- accessors ---
    hh::comm::rank_t rank() const override { return (hh::comm::rank_t)rank_; }
    std::uint32_t nbProcesses() const override { return (std::uint32_t)nbProcesses_; }

    // --- channels ---
    hh::comm::channel_t newChannel() override {
        std::lock_guard<std::mutex> lk(mutex());
        auto channel = (hh::comm::channel_t)comms_.size();
        comms_.push_back(MPI_Comm{});
        MPI_Comm_split(MPI_COMM_WORLD, (int)channel, rank_, &comms_.back());
        return channel;
    }

private:
    int rank_ = -1;
    int nbProcesses_ = -1;
    hh::comm::RequestPool<MPIRequest> requestPool_ = {};
    std::vector<MPI_Comm> comms_ = {};
};

#endif // FDS_COMM_SERVICE_H
