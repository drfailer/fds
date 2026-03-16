#ifndef PRESSURE_ITERATION_DATA_H
#define PRESSURE_ITERATION_DATA_H

#include <memory>
#include <vector>
#include "mesh_data.h"

/// Data token flowing through the pressure iteration sub-graph.
/// Carries mesh tokens and iteration state through the convergence loop.
struct PressureIterData {
    std::vector<std::shared_ptr<MeshData>> meshes;
    double t;
    double dt;

    PressureIterData(const std::vector<std::shared_ptr<MeshData>>& m,
                     double t_, double dt_)
        : meshes(m), t(t_), dt(dt_) {}

    int nm_count() const { return static_cast<int>(meshes.size()); }
};

#endif // PRESSURE_ITERATION_DATA_H
