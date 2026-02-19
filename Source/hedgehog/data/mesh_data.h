#ifndef MESH_DATA_H
#define MESH_DATA_H

#include <ostream>

/// Token type flowing through the Hedgehog dataflow graph.
/// Each token represents one mesh at a given point in the time-stepping pipeline.
struct MeshData {
    int nm;       ///< Mesh index (1-based, Fortran convention)
    double t;     ///< Current simulation time
    double dt;    ///< Current time step
    int phase;    ///< 0 = predictor, 1 = corrector

    MeshData() : nm(0), t(0.0), dt(0.0), phase(0) {}
    MeshData(int nm_, double t_, double dt_, int phase_)
        : nm(nm_), t(t_), dt(dt_), phase(phase_) {}

    friend std::ostream &operator<<(std::ostream &os, const MeshData &md) {
        os << "MeshData{nm=" << md.nm << ", t=" << md.t
           << ", dt=" << md.dt << ", phase=" << md.phase << "}";
        return os;
    }
};

#endif // MESH_DATA_H
