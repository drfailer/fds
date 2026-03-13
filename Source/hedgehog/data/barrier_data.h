#ifndef BARRIER_DATA_H
#define BARRIER_DATA_H

#include <memory>
#include <ostream>
#include <vector>
#include "mesh_data.h"

/// Data type emitted by CollectorState after all mesh tokens arrive at a barrier.
/// Wraps the collected MeshData tokens so that a single-threaded barrier task
/// can perform the global computation and then scatter the tokens back out.
struct BarrierData {
    std::vector<std::shared_ptr<MeshData>> meshes;

    // Fields used by TimestepTask to communicate results to TimestepLoopState
    bool done = false;           ///< True when simulation should terminate
    double newDt = 0.0;          ///< CFL-adjusted DT for next time step
    int newIcyc = 0;             ///< ICYC value for the next time step

    double t() const { return meshes.empty() ? 0.0 : meshes[0]->t; }
    double dt() const { return meshes.empty() ? 0.0 : meshes[0]->dt; }
    int phase() const { return meshes.empty() ? 0 : meshes[0]->phase; }

    friend std::ostream &operator<<(std::ostream &os, const BarrierData &bd) {
        os << "BarrierData{nmeshes=" << bd.meshes.size()
           << ", t=" << bd.t() << ", dt=" << bd.dt() << "}";
        return os;
    }
};

#endif // BARRIER_DATA_H
