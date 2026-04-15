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
    VelErrorPhase         ///< Routing to VelocityErrorTask (single-process)
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
    double dt_bc;      ///< Boundary condition time step (set by WallBC orchestrator barrier)
    int call_ht_1d;    ///< Flag for 1-D heat transfer (0=false, 1=true)
    int exchangeCode;   ///< Exchange operation code (5=flux, 3/6=velocity, 1/4=species)
    int exchangeRound;  ///< Exchange round index (for double-buffered state selection)

    MeshData() : nm(0), t(0.0), dt(0.0), phase(0), firstPass(true), dt_bc(0.0), call_ht_1d(0), exchangeCode(5), exchangeRound(0) {}
    MeshData(int nm_, double t_, double dt_, int phase_)
        : nm(nm_), t(t_), dt(dt_), phase(phase_), firstPass(true), dt_bc(0.0), call_ht_1d(0), exchangeCode(5), exchangeRound(0) {}

    /// Converting constructor: copy fields from a MeshData with a different state tag.
    template<MeshState From>
    explicit MeshData(const MeshData<From>& o)
        : nm(o.nm), t(o.t), dt(o.dt), phase(o.phase), firstPass(o.firstPass),
          dt_bc(o.dt_bc), call_ht_1d(o.call_ht_1d),
          exchangeCode(o.exchangeCode), exchangeRound(o.exchangeRound) {}

    /// Create a shared_ptr<MeshData<To>> with copied fields.
    template<MeshState To>
    std::shared_ptr<MeshData<To>> retag() const {
        return std::make_shared<MeshData<To>>(*this);
    }

    friend std::ostream &operator<<(std::ostream &os, const MeshData &md) {
        os << "MeshData{nm=" << md.nm << ", t=" << md.t
           << ", dt=" << md.dt << ", phase=" << md.phase
           << ", firstPass=" << md.firstPass << "}";
        return os;
    }
};

#endif // MESH_DATA_H
