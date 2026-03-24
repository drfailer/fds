#ifndef FLUX_EXCHANGE_STATE_H
#define FLUX_EXCHANGE_STATE_H

#include <hedgehog/hedgehog.h>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <unordered_map>
#include <vector>
#include "../data/mesh_data.h"
#include "../data/flux_exchange_data.h"
#include "../fds_fortran_interface.h"

// ---------------------------------------------------------------------------
// FluxPackState: receives MeshData, emits FluxExchangeData + MeshData
//
// When MeshData(NM) arrives:
//   1. For each neighbor NOM of NM with NIC_S > 0:
//      a. Call fds_flux_copy_neighbor(NM, NOM) — direct Fortran copy
//      b. Emit FluxExchangeData(sourceNM=NM, destNM=NOM)
//   2. Emit MeshData(NM) directly to FluxCollectorState
//
// This is a state (single-threaded) because mesh arrival order is
// non-deterministic and we need to process each mesh immediately.
// ---------------------------------------------------------------------------

class FluxPackState
    : public hh::AbstractState<1, MeshData, FluxExchangeData, MeshData> {
public:
    explicit FluxPackState(int nmeshes, int nmOffset)
        : nmeshes_(nmeshes), nmOffset_(nmOffset) {
        // Pre-query neighbor topology for each mesh
        for (int i = 0; i < nmeshes; ++i) {
            int nm = nmOffset + i;
            int nNeighbors = fds_flux_get_neighbor_count(nm);
            std::vector<int> neighbors;
            for (int j = 1; j <= nNeighbors; ++j) {
                int nom = fds_flux_get_neighbor_mesh(nm, j);
                if (fds_flux_has_send_cells(nm, nom)) {
                    neighbors.push_back(nom);
                }
            }
            sendNeighbors_[nm] = std::move(neighbors);
        }
    }

    void execute(std::shared_ptr<MeshData> data) override {
        auto t0 = std::chrono::steady_clock::now();
        int nm = data->nm;

        // For each neighbor NOM that NM sends flux data to:
        for (int nom : sendNeighbors_[nm]) {
            // Do the direct Fortran array copy (same-process)
            fds_flux_copy_neighbor(nm, nom);

            // Emit completion token — tells FluxCollectorState that
            // destNM=NOM has received flux data from sourceNM=NM
            this->addResult(std::make_shared<FluxExchangeData>(nm, nom));
        }

        auto t1 = std::chrono::steady_clock::now();
        totalTime_ += std::chrono::duration<double>(t1 - t0).count();
        ++invocations_;

        // Forward MeshData to collector (separate type path)
        this->addResult(data);
    }

    [[nodiscard]] std::string info() const {
        std::ostringstream oss;
        oss << "FLUX_PACK\\n"
            << std::fixed << std::setprecision(3) << totalTime_ << "s"
            << " / " << invocations_ << " calls";
        if (invocations_ > 0)
            oss << " / avg " << std::setprecision(3)
                << (totalTime_ * 1000.0 / invocations_) << "ms";
        return oss.str();
    }

private:
    int nmeshes_, nmOffset_;
    std::unordered_map<int, std::vector<int>> sendNeighbors_;
    double totalTime_ = 0.0;
    int invocations_ = 0;
};

class FluxPackStateManager
    : public hh::StateManager<1, MeshData, FluxExchangeData, MeshData> {
public:
    FluxPackStateManager(std::shared_ptr<FluxPackState> const& state,
                         std::string const& name)
        : hh::StateManager<1, MeshData, FluxExchangeData, MeshData>(state, name) {}

    [[nodiscard]] std::string extraPrintingInformation() const override {
        this->state()->lock();
        auto ret = std::dynamic_pointer_cast<FluxPackState>(this->state())->info();
        this->state()->unlock();
        return ret;
    }
};

// ---------------------------------------------------------------------------
// FluxCollectorState: receives FluxExchangeData + MeshData, emits MeshData
//
// Tracks two things per mesh NM:
//   1. Has MeshData(NM) arrived? (forwarded from FluxPackState)
//   2. How many neighbor FluxExchangeData(dest=NM) have arrived?
//
// When both conditions are met (MeshData present AND all expected neighbors
// have reported), emit MeshData(NM) downstream to the kernel.
//
// This implements the per-neighbor barrier: each mesh waits only for its
// own neighbors, not all meshes.
// ---------------------------------------------------------------------------

class FluxCollectorState
    : public hh::AbstractState<2, FluxExchangeData, MeshData, MeshData> {
public:
    explicit FluxCollectorState(int nmeshes, int nmOffset)
        : nmeshes_(nmeshes), nmOffset_(nmOffset) {
        // Pre-query expected receive count per mesh
        for (int i = 0; i < nmeshes; ++i) {
            int nm = nmOffset + i;
            expectedCount_[nm] = fds_flux_recv_count(nm);
            receivedCount_[nm] = 0;
        }
    }

    void execute(std::shared_ptr<FluxExchangeData> data) override {
        int destNM = data->destNM;
        receivedCount_[destNM]++;
        tryEmit(destNM);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        meshData_[data->nm] = data;
        tryEmit(data->nm);
    }

    [[nodiscard]] std::string info() const {
        std::ostringstream oss;
        oss << "FLUX_COLLECT\\n"
            << emitted_ << " meshes emitted";
        return oss.str();
    }

private:
    void tryEmit(int nm) {
        auto it = meshData_.find(nm);
        if (it == meshData_.end()) return;
        if (receivedCount_[nm] < expectedCount_[nm]) return;

        // All neighbors have reported and MeshData is available
        this->addResult(it->second);
        meshData_.erase(it);
        receivedCount_[nm] = 0;
        ++emitted_;
    }

    int nmeshes_, nmOffset_;
    std::unordered_map<int, int> expectedCount_;
    std::unordered_map<int, int> receivedCount_;
    std::unordered_map<int, std::shared_ptr<MeshData>> meshData_;
    int emitted_ = 0;
};

class FluxCollectorStateManager
    : public hh::StateManager<2, FluxExchangeData, MeshData, MeshData> {
public:
    FluxCollectorStateManager(std::shared_ptr<FluxCollectorState> const& state,
                              std::string const& name)
        : hh::StateManager<2, FluxExchangeData, MeshData, MeshData>(state, name) {}

    [[nodiscard]] std::string extraPrintingInformation() const override {
        this->state()->lock();
        auto ret = std::dynamic_pointer_cast<FluxCollectorState>(this->state())->info();
        this->state()->unlock();
        return ret;
    }
};

#endif // FLUX_EXCHANGE_STATE_H
