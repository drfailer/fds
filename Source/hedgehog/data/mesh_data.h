#ifndef MESH_DATA_H
#define MESH_DATA_H

#include <memory>
#include <ostream>

/// State tag for type-based Hedgehog routing.
/// Each value produces a distinct MeshData<S> type that Hedgehog can route independently.
/// The enum grows as more subgraphs are shared between pipeline phases.
enum class MeshState {
    Default,              ///< Everything outside shared subgraphs
    Init,                 ///< Graph input: initial injection into TimestepState
    PredictorPressure,    ///< Boundary: predictor -> pressure subgraph
    CorrectorPressure,    ///< Boundary: corrector -> pressure subgraph
    Pressure,             ///< Internal pressure pipeline + cycle-back
    SolvePhase,           ///< Routing to PressureSolveKernel (single-process)
    VelErrorPhase,        ///< Routing to VelocityErrorTask (single-process)
    DivExch,              ///< PredDivParallel → DivExchange barrier
    DivP2Pre,             ///< DivExchange barrier → PredDivParallel
    GlobalMat,            ///< PredDivParallel → GlobalMatrix barrier
    DivPart2,             ///< GlobalMatrix barrier → PredDivParallel
    PreSolveExch,         ///< PressureParallel → pre-solve exchange
    PostSolveExch,        ///< PressureParallel → post-solve exchange
    PostParticleOps,      ///< CorrDivSetupCombPartTask → MeshExch7 barrier (after ParticleOps)
    PostVelCorr,          ///< CorrFinalKernelTask → CorrFinalOrch (after velocity correction)
    PostCorrStep1,        ///< CorrStep1 phase → MeshExch4 barrier → DivSetupCombPart phase
    PostHvac,             ///< HvacCalc barrier → CorrDivSetupCombPartTask WallBC phase
    PostWallBC            ///< CorrDivSetupCombPartTask WallBC phase → Fork2 (radiation || DivP1)
};

/// Token type flowing through the Hedgehog dataflow graph.
/// Each token represents one mesh at a given point in the time-stepping pipeline.
/// The template parameter S is a compile-time tag for type-based routing;
/// all instantiations share the same layout and fields.
template<MeshState S = MeshState::Default>
struct MeshData {
    int nm;            ///< Mesh index (1-based, Fortran convention)
    double t;          ///< Current simulation time
    double dt;         ///< Current time step
    int phase;         ///< 0 = predictor, 1 = corrector
    bool firstPass;    ///< True on first pass through CHANGE_TIME_STEP_LOOP, false on CFL retry
    double dt_bc;      ///< Boundary condition time step (computed per-mesh in WallBCKernelTask)
    int call_ht_1d;    ///< Flag for 1-D heat transfer (0=false, 1=true)
    int wall_counter;  ///< Per-mesh copy of global WALL_COUNTER (incremented each corrector step)
    int exchangeCode;   ///< Exchange operation code (5=flux, 3/6=velocity, 1/4=species)
    int exchangeRound;  ///< Exchange round index (for double-buffered state selection)
    MeshData() : nm(0), t(0.0), dt(0.0), phase(0), firstPass(true), dt_bc(0.0), call_ht_1d(0), wall_counter(0), exchangeCode(5), exchangeRound(0) {}
    MeshData(int nm_, double t_, double dt_, int phase_)
        : nm(nm_), t(t_), dt(dt_), phase(phase_), firstPass(true), dt_bc(0.0), call_ht_1d(0), wall_counter(0), exchangeCode(5), exchangeRound(0) {}

    friend std::ostream &operator<<(std::ostream &os, const MeshData &md) {
        os << "MeshData{nm=" << md.nm << ", t=" << md.t
           << ", dt=" << md.dt << ", phase=" << md.phase
           << ", firstPass=" << md.firstPass << "}";
        return os;
    }
};

/// Zero-cost retag: reinterpret the same MeshData object as a different state tag.
/// Safe because all MeshData<S> instantiations share identical layout.
template<MeshState To, MeshState From>
std::shared_ptr<MeshData<To>> retag(std::shared_ptr<MeshData<From>> p) {
    return std::reinterpret_pointer_cast<MeshData<To>>(std::move(p));
}

#endif // MESH_DATA_H
