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

## Architecture: Push-Then-Gate

The global `MESH_EXCHANGE(CODE)` barrier is decomposed into two nodes:

```
[ExchangePushTask]    parallel task, copies source data to targets' OMESHes
        |
[ExchangeGateState]   single-threaded state, pure dependency tracker
```

### ExchangePushTask (parallel, numThreads = kernelThreads)

When mesh NM arrives, the push task copies NM's boundary data to all same-rank
target meshes.  For CODE 5 (momentum fluxes):

```cpp
void execute(std::shared_ptr<MeshData> md) override {
    for (int target : depGraph_->sendTargets(md->nm)) {
        if (fds_mesh_process(target) == myRank_) {
            fds_flux_copy_neighbor_ts(md->nm, target);
        }
    }
    this->addResult(md);   // signal "push done"
}
```

`fds_flux_copy_neighbor_ts(NM, NOM)` reads from `MESHES(NM)` and writes to
`MESHES(NOM)%OMESH(NM)`.  This is a **push** operation: NM pushes its data
out to NOM's ghost buffer.

Thread safety of parallel pushes:
- Different source meshes write to different OMESH entries of the target
  (`OMESH(A)` vs `OMESH(C)`), so there are no write conflicts.
- Reads from source meshes are concurrent-safe (read-only before any solve).

### ExchangeGateState (single-threaded, pure dependency tracker)

Receives "push done" signals from the upstream push task.  Each signal means
mesh NM has finished copying its data to all targets.  The gate tracks which
meshes have pushed and emits a mesh downstream only when **all of its
receive-dependencies** have also pushed:

```cpp
void execute(std::shared_ptr<MeshData> data) override {
    int nm = data->nm;
    pendingMeshes_[nm] = data;
    satisfied_.set(nm - 1);

    // Check all local meshes for newly-satisfied dependencies
    for (int destNM = lower_; destNM <= upper_; ++destNM) {
        if (emitted_[destNM]) continue;
        if (!pendingMeshes_[destNM]) continue;
        if (noDeps_[destNM] ||
            satisfied_.containsAll(depGraph_->sameRankRecvDeps(destNM))) {
            emitted_[destNM] = true;
            this->addResult(std::move(pendingMeshes_[destNM]));
        }
    }
    if (allEmitted()) resetRound();
}
```

A mesh can proceed when:
1. It has been pushed (its data is saved in all targets' OMESHes).
2. All of its receive-dependencies have been pushed (its own OMESH has
   fresh data from all sources).

The gate contains **no I/O** -- all copies and communication happen in the
upstream push task.  This makes the gate reusable across different exchange
codes (5, 3, 1, etc.) by pairing it with different push tasks.

## Why Push-Then-Gate Avoids Data Races

A naive pull model -- where each mesh pulls neighbor data into itself on
arrival -- creates a race condition:

```
1. Pull model processes A: copies B's data into A's OMESH(B).  Emits A.
2. A goes to solve, which modifies A's FVX/FVY/FVZ/H.
3. Pull model processes B: copies A's data into B's OMESH(A).
   But step 2 is modifying A concurrently!  DATA RACE.
```

The push model avoids this because the copy direction is reversed:

```
1. Push processes A: copies A's data TO targets' OMESHes.  Reads A (safe:
   A hasn't entered the solve yet).  Emits A as "push done."
2. Push processes B: copies B's data TO targets' OMESHes (including A's
   OMESH(B)).  Reads B (safe: B hasn't entered the solve yet).  Emits B.
3. Gate receives A and B.  Both have pushed.  Gate releases both.
4. A and B enter the solve.  All OMESH entries are already populated.
   No concurrent reads/writes on source mesh data.
```

Key invariant: **each mesh's push reads from itself before it is released to
the solve**.  No other push task reads from the same mesh, so no solve can
race with a push read.

The gate provides the ordering guarantee: it holds a mesh until all its
sources have pushed (and therefore finished writing to its OMESH entries).
Hedgehog's queue synchronization ensures memory visibility between the push
thread and the solve thread.

## Supporting Infrastructure

### DynBitset (`tool/dyn_bitset.h`)

Dynamic bitset for tracking which meshes have pushed.  O(1) `set()`,
`contains()`, and `containsAll()` (subset check via word-level AND):

```cpp
class DynBitset {
    size_t nbits_;
    std::vector<uint64_t> words_;
public:
    explicit DynBitset(size_t nbits);
    void set(size_t pos);
    bool contains(size_t pos) const;
    bool containsAll(const DynBitset &other) const;  // is other a subset?
    size_t count() const;
    void reset();
};
```

`containsAll()` is the critical operation: it checks whether all bits set in
`other` are also set in `this`.  For N meshes with W = ceil(N/64) words, the
check is O(W) -- typically 1 word for up to 64 meshes.

### MeshDependencyGraph (`tool/mesh_dependency_graph.h`)

Pre-computed dependency graph built once at initialization from Fortran
NIC_R/NIC_S topology queries:

```cpp
class MeshDependencyGraph {
    std::vector<DynBitset> recvDeps_;         // who sends TO mesh nm
    std::vector<DynBitset> sameRankRecvDeps_; // same-rank subset
    std::vector<std::vector<int>> sendTargets_; // who nm sends TO
public:
    MeshDependencyGraph(int lowerMesh, int upperMesh);
    const DynBitset &sameRankRecvDeps(int nm) const;
    const std::vector<int> &sendTargets(int nm) const;
    int totalMeshes() const;
};
```

Construction queries:
- `fds_exchange_recv_dep_count(nm)` / `fds_exchange_recv_dep_mesh(nm, idx)`:
  meshes that send TO nm (`MESHES(NM)%OMESH(NOM)%NIC_R > 0`).
- `fds_exchange_send_dep_count(nm)` / `fds_exchange_send_dep_mesh(nm, idx)`:
  meshes that nm sends TO (`MESHES(NM)%OMESH(NOM)%NIC_S > 0`).
- `fds_mesh_process(nm)`: MPI rank owning mesh nm.

The graph is shared (read-only after construction) between predictor and
corrector subgraphs.

### Mesh Topology Examples

**Abutting meshes** (A shares a face with B):
- `sendTargets(A)` includes B; `sendTargets(B)` includes A (symmetric).
- `recvDeps(A)` includes B; `recvDeps(B)` includes A (symmetric).
- Gate holds both until the other has pushed.

**Embedded meshes** (M3 is coarse, M5 is fine, nested inside M3):
- `sendTargets(M3)` includes M5; `sendTargets(M5)` does NOT include M3.
- `recvDeps(M5)` includes M3; `recvDeps(M3)` does NOT include M5.
- Gate holds M5 until M3 pushes.  M3 can proceed independently.
- This correctly reflects the unidirectional data flow: coarse mesh
  boundary data flows to the fine mesh, not the other way around (for CODE 5).

**Non-neighboring meshes**:
- `sendTargets(A)` does not include C; `recvDeps(A)` does not include C.
- Gate releases A without waiting for C.  This is the key advantage over
  the global barrier.

### Thread-Safe Fortran Routine (`fds_driver.f90`)

The per-neighbor copy routine uses local pointer aliases instead of
module-level M/M2/M3 to enable concurrent calls from different threads:

```fortran
RECURSIVE SUBROUTINE MESH_EXCHANGE_FLUX_NEIGHBOR_TS(NM, NOM)
INTEGER, INTENT(IN) :: NM, NOM
TYPE(MESH_TYPE), POINTER :: ML
TYPE(OMESH_TYPE), POINTER :: OM_SEND, OM_RECV

ML => MESHES(NM)
OM_SEND => ML%OMESH(NOM)
IF (OM_SEND%NIC_S == 0) RETURN

OM_RECV => MESHES(NOM)%OMESH(NM)

OM_RECV%FVX(IMIN:IMAX,...) = ML%FVX(IMIN:IMAX,...)
OM_RECV%FVY(IMIN:IMAX,...) = ML%FVY(IMIN:IMAX,...)
OM_RECV%FVZ(IMIN:IMAX,...) = ML%FVZ(IMIN:IMAX,...)
OM_RECV%H(IMIN:IMAX,...)   = ML%H(IMIN:IMAX,...)   ! or HS
END SUBROUTINE
```

Thread safety: two concurrent calls `(A, B)` and `(C, D)` write to
different OMESH entries of different meshes.  Reads from source meshes are
concurrent-safe.

## Integration: Pressure Iteration Subgraph

The push-then-gate replaces the global `fds_mesh_exchange(5)` barrier in the
first pressure iteration.  The pressure iteration subgraph pipeline:

```
BaroclinicKernel (parallel, N threads)
    |
ExchangePush (parallel, N threads: copy to targets' OMESHes)
    |
ExchangeGate (state: pure dependency tracker)
    |
PressureSolve (parallel, N threads: FFT or ULMAT)
    |
PressureConvergence (barrier: exchange(5) + velocity error + check)
    |--- not converged ---> PressureIterMeshData ---> BaroclinicKernel (cycle)
    |--- converged -------> MeshData ---> subgraph output
```

The BaroclinicKernel computes the baroclinic correction per mesh, then each
mesh's flux data is pushed to its targets.  The gate releases each mesh to
the pressure solve as soon as its dependencies are met.  For non-neighboring
meshes, this can happen before all meshes finish the baroclinic kernel --
overlapping the slow mesh's kernel with the fast mesh's solve.

### Wiring

```cpp
auto depGraph = std::make_shared<MeshDependencyGraph>(lower, upper);

auto baroclinicKernel = std::make_shared<BaroclinicKernelTask>(kernelThreads);
auto exchangePushTask = std::make_shared<ExchangePushTask>(kernelThreads, depGraph);
auto exchangeGateSM   = std::make_shared<ExchangeGateManager>(
    std::make_shared<ExchangeGateState>(depGraph), "ExchangeGate");
auto solveKernel      = std::make_shared<PressureSolveKernelTask>(kernelThreads, presFlag);
auto convergenceSM    = std::make_shared<PressureConvergenceManager>(...);

subgraph->inputs(baroclinicKernel);
subgraph->edges(baroclinicKernel, exchangePushTask);
subgraph->edges(exchangePushTask, exchangeGateSM);
subgraph->edges(exchangeGateSM, solveKernel);
subgraph->edges(solveKernel, convergenceSM);
subgraph->edge<PressureIterMeshData>(convergenceSM, baroclinicKernel);
subgraph->outputs(convergenceSM);
```

The MeshDependencyGraph is built once in `fds_graph.h` and shared by both
predictor and corrector pressure iteration subgraphs.

## Extending to Other Exchange Codes

The ExchangeGateState is code-agnostic -- it tracks dependencies using
bitsets without knowing what data is being exchanged.  To replace a different
`MESH_EXCHANGE(CODE)` barrier:

1. **Create a thread-safe copy routine** for the specific CODE, following the
   `MESH_EXCHANGE_FLUX_NEIGHBOR_TS` pattern (RECURSIVE, local pointers,
   no module-level state).

2. **Create a new push task** that calls the code-specific copy routine:
   ```cpp
   class ExchangePushCodeXTask
       : public hh::AbstractTask<1, MeshData, MeshData> {
       void execute(std::shared_ptr<MeshData> md) override {
           for (int target : depGraph_->sendTargets(md->nm)) {
               if (fds_mesh_process(target) == myRank_) {
                   fds_code_x_copy_neighbor_ts(md->nm, target);
               }
           }
           this->addResult(md);
       }
   };
   ```

3. **Reuse the same ExchangeGateState** with the same MeshDependencyGraph
   (or a code-specific one if the dependency topology differs by CODE).

4. **Wire**: `upstream → PushTask → GateState → downstream`.

For cross-rank exchange (MPI), the push task would additionally pack data
into a send buffer and post `MPI_Isend`.  The gate state would wait for
both same-rank pushes AND `MPI_Irecv` completions before releasing a mesh.

## Files

```
Source/hedgehog/
  tool/dyn_bitset.h                  Dynamic bitset for dependency tracking
  tool/mesh_dependency_graph.h       Pre-computed mesh exchange topology
  task/exchange_push_task.h          Parallel push task (CODE 5 flux copy)
  state/exchange_orchestrator_state.h ExchangeGateState + ExchangeGateManager
  graph/pressure_iteration_subgraph.h Wiring (push + gate in pipeline)
  graph/fds_graph.h                  MeshDependencyGraph construction
  fds_driver.f90                     MESH_EXCHANGE_FLUX_NEIGHBOR_TS
  fds_c_interface.f90                C bindings for dependency queries + copy
  fds_fortran_interface.h            C declarations
```
