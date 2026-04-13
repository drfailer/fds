#ifndef PRESSURE_ITERATION_DATA_H
#define PRESSURE_ITERATION_DATA_H

#include "mesh_data.h"

/// Convenience aliases for pressure iteration MeshData types.
using PressureMeshData      = MeshData<MeshState::Pressure>;
using SolvePhaseMeshData    = MeshData<MeshState::SolvePhase>;
using VelErrorPhaseMeshData = MeshData<MeshState::VelErrorPhase>;
using PredPressureMeshData  = MeshData<MeshState::PredictorPressure>;
using CorrPressureMeshData  = MeshData<MeshState::CorrectorPressure>;

#endif // PRESSURE_ITERATION_DATA_H
