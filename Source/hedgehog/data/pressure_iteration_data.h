#ifndef PRESSURE_ITERATION_DATA_H
#define PRESSURE_ITERATION_DATA_H

#include <memory>
#include <vector>
#include "mesh_data.h"

/// MeshData wrapper for the pressure iteration cycle.
struct PressureIterMeshData {
    std::shared_ptr<MeshData> mesh;
};

/// Wrapper type for routing MeshData to PressureSolveKernel after exchange.
struct SolvePhaseData {
    std::shared_ptr<MeshData> mesh;
};

/// Wrapper type for routing MeshData to VelocityError after exchange.
struct VelErrorPhaseData {
    std::shared_ptr<MeshData> mesh;
};

#endif // PRESSURE_ITERATION_DATA_H
