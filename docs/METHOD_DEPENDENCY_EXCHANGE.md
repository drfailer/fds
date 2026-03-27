# Dependency-Aware Mesh Exchange

## Problem

FDS meshes exchange boundary data via `MESH_EXCHANGE(CODE)`, a global barrier
that synchronizes all meshes before any can proceed.  In the original loop:

```fortran
DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   ! per-mesh work
ENDDO
CALL MESH_EXCHANGE(5)     ! <-- global barrier: waits for ALL meshes
DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   ! per-mesh work that reads OMESH
ENDDO
```

Inside a Hedgehog dataflow graph, this barrier forces all N meshes to complete
one phase before any mesh can start the next.  For N meshes on P threads, the
barrier idle time is `(P - 1) * T_last`, where `T_last` is the time the
slowest mesh takes to arrive.

The goal is to replace the global barrier with per-mesh dependency tracking so
that a mesh can proceed as soon as its specific dependencies are met, without
waiting for unrelated meshes.

## Architecture: Pull-Only Exchange with Dependency Manager

The global `MESH_EXCHANGE(CODE)` barrier is decomposed into a reusable
sub-graph (`MeshExchangeGraph`) containing a cycle between a state machine
and a parallel exchange task:

```
MeshDepsManager (state) <---> FluxExchangeTask (task, parallel)
        |                             ^
  input T tokens                MeshExchangeData (cycle)
        |
  output T tokens (Done meshes)
```

### MeshDependenciesManagerState (single-threaded state, dependency tracker)

Manages a 4-state machine per mesh:

```
NotArrived -> Wait -> Processing -> Processed -> Done
```

**NotArrived -> Wait**: When a mesh token arrives from the upstream kernel,
the mesh enters Wait and decrements `unarrivedNeighborCount` for all its
same-rank neighbors (notifying them that one more dependency has arrived).

**Wait -> Processing**: When `unarrivedNeighborCount[nm] == 0` (all neighbors
have arrived and are therefore safe to read from), the state emits a
`MeshExchangeData` token containing the mesh and its full neighbor list.
This dispatches the mesh to the exchange task.

**Processing -> Processed**: When the exchange task returns the
`MeshExchangeData` (cycle edge), the mesh enters Processed and decrements
`unprocessedNeighborCount` for all its neighbors.

**Processed -> Done**: When `unprocessedNeighborCount[nm] == 0` (all neighbors
have finished their exchange and are therefore not reading from this mesh
anymore), the original `MeshData` token is emitted downstream.

The Done gate is critical: it prevents a mesh from entering the solve (which
modifies source arrays like FVX/FVY/FVZ/H) while a neighbor is still pulling
data from it in the exchange task.

Counter-based dependency tracking avoids full-mesh scans:
- `unarrivedNeighborCount[nm]`: initialized to neighbor count, decremented on
  each neighbor's arrival.  Guards Wait -> Processing.
- `unprocessedNeighborCount[nm]`: initialized to neighbor count, decremented on
  each neighbor's exchange return.  Guards Processed -> Done.

### FluxExchangeTask (parallel task, numThreads from budget)

Performs pull-only copies: each mesh pulls data **from** all its same-rank
neighbors **into** its own OMESH buffers.  For CODE 5 (momentum fluxes):

```cpp
void execute(std::shared_ptr<MeshExchangeData> data) override {
    int nm = data->mesh->nm;
    for (int nom : data->neighbors) {
        fds_flux_copy_neighbor_ts(nom, nm);  // pull NOM's data into NM's OMESH
    }
    this->addResult(data);
}
```

`fds_flux_copy_neighbor_ts(NOM, NM)` reads from `MESHES(NOM)` and writes to
`MESHES(NM)%OMESH(NOM)`.

### MeshExchangeGraph (reusable sub-graph)

Encapsulates the dependency manager + exchange task cycle as a typed sub-graph
that can be instantiated multiple times in a pipeline:

```cpp
template <typename T>
class MeshExchangeGraph : public hh::Graph<2, T, TerminationData, T> {
    MeshExchangeGraph(depGraph, exchangeTask, name) {
        // Entry: T -> depManager, TerminationData -> depManager
        // Cycle: depManager -> exchangeTask -> depManager
        // Exit: depManager emits T (Done meshes)
    }
};
```

## Why Pull-Only Avoids Data Races

Each mesh writes only to its **own** OMESH entries.  Different meshes writing
concurrently never conflict because they target different OMESH arrays:

```
Mesh A pulls from B: writes to MESHES(A)%OMESH(B)
Mesh B pulls from A: writes to MESHES(B)%OMESH(A)
```

These are entirely separate memory locations.  No synchronization is needed
between concurrent pull operations.

The dependency manager provides two ordering guarantees:

1. **Wait -> Processing gate**: A mesh only starts pulling from neighbors
   after all neighbors have arrived (are at least in Wait state).  This
   ensures the exchange reads consistent source data — no neighbor is still
   being modified by a previous-phase kernel.

2. **Processed -> Done gate**: A mesh only proceeds to the next kernel after
   all its neighbors have finished pulling.  This ensures no neighbor is
   still reading the mesh's source data while the next kernel modifies it.

```
Timeline for meshes A (fast) and B (slow):

1. A arrives, enters Wait.  B has not arrived yet.
   unarrivedNeighborCount[A] = 1 (waiting for B).

2. B arrives, enters Wait.  Decrements A's counter to 0.
   Both A and B transition Wait -> Processing.

3. A's exchange task pulls from B.  B's exchange task pulls from A.
   (Concurrent: write to different OMESHs.)

4. A returns, enters Processed.  Decrements B's unprocessedNeighborCount.
   A cannot enter Done yet: B is still in Processing (reading from A).

5. B returns, enters Processed.  Decrements A's unprocessedNeighborCount.
   Both counters reach 0 -> both transition to Done.

6. A and B are emitted downstream.  Safe to start the solve.
```

## Thread Budget

Exchange task threads are allocated from the global `ThreadBudget`:

```cpp
auto preSolveExchange = std::make_shared<MeshExchangeGraph<MeshData>>(
    depGraph, std::make_shared<FluxExchangeTask>(budget.fluxExchange),
    "PreSolveExchange");
```

The budget distributes threads proportionally by weight across pipeline stages.
Exchange tasks are LIGHT weight (short per-element copy) relative to kernel
tasks, so they typically receive fewer threads than the pressure solve.

Within the pressure iteration pipeline, threads are distributed:

```
Baroclinic(1) -> Exchange(1) -> Solve(4) -> Exchange(1) -> VelError(2)
```

Weights sum to 8.  On a machine with 40 hardware threads and 16 meshes
(cap=16), this yields: baroclinic=2, exchange=2, solve=8, velError=4.

## Supporting Infrastructure

### MeshExchangeData (`data/mesh_exchange_data.h`)

Token flowing through the exchange cycle:

```cpp
struct MeshExchangeData {
    std::shared_ptr<MeshData> mesh;   // The mesh to exchange
    std::vector<int> neighbors;       // Same-rank neighbors to pull from (1-based)
};
```

### MeshDependencyGraph (`tool/mesh_dependency_graph.h`)

Pre-computed dependency graph built once at initialization from Fortran
NIC_R/NIC_S topology queries:

```cpp
class MeshDependencyGraph {
    std::vector<DynBitset> recvDeps_;               // who sends TO mesh nm
    std::vector<DynBitset> sameRankRecvDeps_;       // same-rank subset
    std::vector<std::vector<int>> sendTargets_;     // who nm sends TO
    std::vector<DynBitset> sameRankNeighbors_;      // union(recv, send) same-rank
    std::vector<std::vector<int>> sameRankNeighborsList_;  // sorted list form
public:
    MeshDependencyGraph(int lowerMesh, int upperMesh);
    const std::vector<int> &sameRankNeighborsList(int nm) const;
};
```

The `sameRankNeighborsList` is the union of receive dependencies and send
targets, filtered to the same MPI rank.  This is what the dependency manager
uses: a mesh must wait for ALL its neighbors (both those it reads from and
those that read from it) because the Done gate needs the symmetry.

Construction queries:
- `fds_exchange_recv_dep_count(nm)` / `fds_exchange_recv_dep_mesh(nm, idx)`:
  meshes that send TO nm (`MESHES(NM)%OMESH(NOM)%NIC_R > 0`).
- `fds_exchange_send_dep_count(nm)` / `fds_exchange_send_dep_mesh(nm, idx)`:
  meshes that nm sends TO (`MESHES(NM)%OMESH(NOM)%NIC_S > 0`).
- `fds_mesh_process(nm)`: MPI rank owning mesh nm.

The graph is shared (read-only after construction) between predictor and
corrector pressure iteration subgraphs.

### DynBitset (`tool/dyn_bitset.h`)

Dynamic bitset for tracking which meshes have arrived/processed.  O(1) `set()`,
`contains()`, and `containsAll()` (subset check via word-level AND):

```cpp
class DynBitset {
    std::vector<uint64_t> words_;
public:
    explicit DynBitset(size_t nbits);
    void set(size_t pos);
    bool contains(size_t pos) const;
    bool containsAll(const DynBitset &other) const;
    void reset();
};
```

### Mesh Topology Examples

**Abutting meshes** (A shares a face with B):
- `sameRankNeighborsList(A)` includes B; `sameRankNeighborsList(B)` includes A.
- Both must arrive before either starts exchanging.
- Both must finish exchanging before either proceeds to solve.

**Embedded meshes** (M3 coarse, M5 fine, nested inside M3):
- `sameRankNeighborsList(M3)` includes M5; `sameRankNeighborsList(M5)` includes M3.
- Both wait for each other (symmetric neighbors), ensuring the Done gate
  works correctly in both directions.

**Non-neighboring meshes**:
- `sameRankNeighborsList(A)` does not include C.
- The dependency manager releases A without waiting for C.  This is the key
  advantage over the global barrier.

### Thread-Safe Fortran Routine (`fds_driver.f90`)

The per-neighbor copy routine uses local pointer aliases instead of
module-level M/M2/M3 to enable concurrent calls from different threads:

```fortran
RECURSIVE SUBROUTINE MESH_EXCHANGE_FLUX_NEIGHBOR_TS(NOM, NM)
INTEGER, INTENT(IN) :: NOM, NM
TYPE(MESH_TYPE), POINTER :: ML
TYPE(OMESH_TYPE), POINTER :: OM_SEND, OM_RECV

ML => MESHES(NOM)
OM_SEND => ML%OMESH(NM)
IF (OM_SEND%NIC_S == 0) RETURN

OM_RECV => MESHES(NM)%OMESH(NOM)

OM_RECV%FVX(IMIN:IMAX,...) = ML%FVX(IMIN:IMAX,...)
OM_RECV%FVY(IMIN:IMAX,...) = ML%FVY(IMIN:IMAX,...)
OM_RECV%FVZ(IMIN:IMAX,...) = ML%FVZ(IMIN:IMAX,...)
OM_RECV%H(IMIN:IMAX,...)   = ML%H(IMIN:IMAX,...)   ! or HS
END SUBROUTINE
```

Thread safety: calls `(A pulling from B)` and `(C pulling from D)` write to
different OMESH entries of different meshes.  Reads from source meshes are
concurrent-safe (source data is immutable until the Done gate releases the mesh).

## Integration: Pressure Iteration Subgraph

The pull-only exchange replaces the global `fds_mesh_exchange(5)` in the
pressure iteration.  Two exchange instances are used — one before and one
after the pressure solve:

```
BaroclinicKernel (parallel)
    |
PreSolveExchange (MeshExchangeGraph: pull-only parallel exchange)
    |
PressureSolve (parallel: FFT or ULMAT)
    |
PostSolveExchange (MeshExchangeGraph: pull-only parallel exchange)
    |
VelocityErrorTask (parallel)
    |
PressureConvergence (barrier: convergence check only)
    |--- not converged ---> PressureIterMeshData ---> BaroclinicKernel (cycle)
    |--- converged -------> MeshData ---> subgraph output
```

### Wiring

```cpp
auto depGraph = std::make_shared<MeshDependencyGraph>(lower, upper);

auto baroclinicKernel = std::make_shared<BaroclinicKernelTask>(budget.baroclinic);
auto solveKernel = std::make_shared<PressureSolveKernelTask>(budget.pressureSolve, presFlag);
auto velErrorTask = std::make_shared<VelocityErrorTask>(budget.velError);

auto preSolveExchange = std::make_shared<MeshExchangeGraph<MeshData>>(
    depGraph, std::make_shared<FluxExchangeTask>(budget.fluxExchange),
    "PreSolveExchange");
auto postSolveExchange = std::make_shared<MeshExchangeGraph<MeshData>>(
    depGraph, std::make_shared<FluxExchangeTask>(budget.fluxExchange),
    "PostSolveExchange");

auto convergenceSM = std::make_shared<PressureConvergenceManager>(...);

subgraph->input<MeshData>(baroclinicKernel);
subgraph->edges(baroclinicKernel, preSolveExchange);
subgraph->edges(preSolveExchange, solveKernel);
subgraph->edges(solveKernel, postSolveExchange);
subgraph->edges(postSolveExchange, velErrorTask);
subgraph->edges(velErrorTask, convergenceSM);
subgraph->edge<PressureIterMeshData>(convergenceSM, baroclinicKernel);
subgraph->outputs(convergenceSM);
```

The MeshDependencyGraph is built once in `fds_graph.h` and shared by both
predictor and corrector pressure iteration subgraphs.

### Termination

Each `MeshExchangeGraph` instance contains an internal cycle (state <-> task).
To terminate cleanly when the simulation ends, `TerminationData` is routed
from the parent graph input to each exchange sub-graph:

```cpp
MeshExchangeGraph<MeshData>::wireTermination(subgraph, preSolveExchange);
MeshExchangeGraph<MeshData>::wireTermination(subgraph, postSolveExchange);
```

When `TerminationData` arrives, the dependency manager sets `done_=true` and
its `canTerminate()` override returns `true`, allowing the cycle to shut down.

## Extending to Other Exchange Codes

The `MeshExchangeGraph` is exchange-code-agnostic.  The dependency manager
tracks arrival/completion states without knowing what data is copied.  To
replace a different `MESH_EXCHANGE(CODE)` barrier:

1. **Create a thread-safe copy routine** for the specific CODE, following the
   `MESH_EXCHANGE_FLUX_NEIGHBOR_TS` pattern (RECURSIVE, local pointers,
   no module-level state).

2. **Create a new exchange task** that calls the code-specific copy routine:
   ```cpp
   class CodeXExchangeTask
       : public hh::AbstractTask<1, MeshExchangeData, MeshExchangeData> {
       void execute(std::shared_ptr<MeshExchangeData> data) override {
           int nm = data->mesh->nm;
           for (int nom : data->neighbors) {
               fds_code_x_copy_neighbor_ts(nom, nm);
           }
           this->addResult(data);
       }
   };
   ```

3. **Instantiate a MeshExchangeGraph** with the new task:
   ```cpp
   auto exchange = std::make_shared<MeshExchangeGraph<MeshData>>(
       depGraph, std::make_shared<CodeXExchangeTask>(budget.someField),
       "CodeXExchange");
   ```

4. **Wire** in the parent graph and route `TerminationData`.

For cross-rank exchange (MPI), the exchange task would additionally pack data
into a send buffer and post `MPI_Isend`.  The dependency manager would need
additional state to track `MPI_Irecv` completions before transitioning to Done.

## Files

```
Source/hedgehog/
  data/mesh_exchange_data.h               Token: mesh + neighbor list
  tool/dyn_bitset.h                       Dynamic bitset for dependency tracking
  tool/mesh_dependency_graph.h            Pre-computed mesh exchange topology
  tool/thread_budget.h                    Hardware-aware thread allocation
  task/mesh_exchange_task.h               FluxExchangeTask (parallel pull-only)
  task/velocity_error_task.h              VelocityErrorTask (parallel)
  state/mesh_dependencies_manager_state.h MeshDepsManager state + StateManager
  graph/mesh_exchange_graph.h             Reusable exchange sub-graph template
  graph/pressure_iteration_subgraph.h     Pipeline wiring (2 exchange instances)
  graph/fds_graph.h                       MeshDependencyGraph construction
  fds_driver.f90                          MESH_EXCHANGE_FLUX_NEIGHBOR_TS
  fds_c_interface.f90                     C bindings for dependency queries + copy
  fds_fortran_interface.h                 C declarations
```
