# Pattern B Complex Routine Parallelization

## Overview

This methodology covers parallelizing large, complex routines with cross-mesh dependencies using Hedgehog's Pattern B architecture. Pattern B uses **sequential preprocessing and finalization** with **parallel kernel execution** in the middle.

**When to use Pattern B Complex**:
- Routine has 150+ lines of orchestration code
- Contains cross-mesh dependencies (OMESH reads/writes)
- Multiple phases with different dependency patterns
- Callees need thread-safe conversion
- Only a portion of the routine is parallelizable

**Example**: WALL_BC (239 lines) → WallBC sub-graph (~90% parallelizable)

## Architecture: Three-Phase Pattern

```
┌─────────────────────────────────────────────────────────┐
│  ORCHESTRATOR (sequential, collects N meshes)           │
│  1. Compute global parameters from shared state         │
│  2. Run sequential preprocessing (OMESH reads)           │
│  3. Dispatch N parallel work tokens                      │
└─────────────────────────────────────────────────────────┘
                         ↓
┌─────────────────────────────────────────────────────────┐
│  KERNEL TASK (parallel, N threads)                      │
│  - Process ~80-90% of work without cross-mesh access    │
│  - Thread-safe kernels only                              │
└─────────────────────────────────────────────────────────┘
                         ↓
┌─────────────────────────────────────────────────────────┐
│  COLLECTOR (sequential, gathers N results)              │
│  1. Sort by mesh index (deterministic ordering)         │
│  2. Run sequential finalization (OMESH writes)           │
│  3. Emit N MeshData tokens                               │
└─────────────────────────────────────────────────────────┘
```

### Phase 1: Sequential Preprocessing (Orchestrator)

**Purpose**: Handle operations that require cross-mesh coordination or global state

**Common operations**:
- Compute global parameters (time steps, flags, counters)
- Read from neighboring meshes (OMESH, INTERPOLATED_BOUNDARY)
- Setup shared data structures
- Initialize per-mesh state before parallel execution

**Example from WallBC**:
```cpp
// Compute global parameters once per time step
double dt_bc = fds_compute_wall_bc_dt_bc(collected_[0]->t);
int call_ht_1d = fds_check_call_ht_1d();

// Update global state if needed
if (call_ht_1d) {
    fds_update_bc_clock(collected_[0]->t);
}

// Sequential preprocessing for each mesh
for (auto &md : collected_) {
    fds_wall_bc_preprocessing(md->nm, md->t, dt_bc, call_ht_1d);
}
```

### Phase 2: Parallel Kernel Execution

**Purpose**: Process the bulk of work without cross-mesh dependencies

**Characteristics**:
- Operates on local mesh data only (TYPE(MESH_TYPE) argument)
- No OMESH access
- No global state modification
- Thread-safe (RECURSIVE keyword in Fortran)
- Typically 80-90% of total work

**Example from WallBC**:
```fortran
RECURSIVE SUBROUTINE WALL_BC_PROCESS_CELLS_KERNEL(M, NM, PREDICTOR_FLAG, T, DT, DT_BC, CALL_HT_1D)
  TYPE(MESH_TYPE), POINTER :: M

  ! Process wall cells WITHOUT cross-mesh flags
  DO IW = 1, N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS
    WC => M%WALL(IW)

    ! Skip cells requiring sequential processing
    IF (WC%HAS_INTERPOLATED_BC) CYCLE
    IF (WC%HAS_BACK_MESH) CYCLE

    ! Thread-safe processing
    CALL SURFACE_HEAT_TRANSFER(NM, PREDICTOR_FLAG, ...)
    CALL CALCULATE_ZZ_F(NM, PREDICTOR_FLAG, ...)
  END DO
END SUBROUTINE
```

### Phase 3: Sequential Finalization (Collector)

**Purpose**: Handle operations requiring cross-mesh writes or mesh ordering

**Common operations**:
- Sort results by mesh index (deterministic output)
- Write to neighboring meshes (OMESH updates)
- Handle cross-mesh coupling (BACK_MESH, particle transfer)
- Finalize shared data structures

**Example from WallBC**:
```cpp
// Sort results deterministically
std::sort(results_.begin(), results_.end(),
          [](const auto &a, const auto &b) { return a->nm < b->nm; });

// Sequential finalization for each mesh
for (auto &w : results_) {
    fds_wall_bc_finalize(w->nm, w->t, w->dt_bc, w->call_ht_1d);
}
```

## Step-by-Step Implementation

### Step 1: Analyze the Original Routine

Identify the different phases and their dependencies:

```fortran
SUBROUTINE WALL_BC(T, DT, NM)
  ! Phase 1: Cross-mesh reads (OMESH access)
  IF (N_EXTERNAL_WALL_CELLS > 0) CALL ASSIGN_GHOST_VALUE(...)  ! OMESH

  ! Common setup (can be in preprocessing)
  DO IW = 1, N_WALL_CELLS
    CALL NEAR_SURFACE_GAS_VARIABLES(...)  ! Local only
    IF (CALL_HT_1D) THEN
      HEAT_COEF = HEAT_TRANSFER_COEFFICIENT(...)  ! Local only
    ENDIF
  END DO

  ! Phase 2: Main cell processing (90% of cells, local only)
  DO IW = 1, N_WALL_CELLS
    IF (.NOT. HAS_INTERPOLATED_BC .AND. .NOT. HAS_BACK_MESH) THEN
      CALL SURFACE_HEAT_TRANSFER(...)  ! Local only
      CALL CALCULATE_ZZ_F(...)  ! Local only
    ENDIF
  END DO

  ! Phase 3: Cross-mesh writes (OMESH access)
  DO IW = 1, N_WALL_CELLS
    IF (HAS_BACK_MESH) THEN
      CALL SOLID_HEAT_TRANSFER(..., BACK_MESH=...)  ! Cross-mesh coupling
    ENDIF
  END DO

  IF (CORRECTOR) CALL DEPOSIT_PARTICLE_MASS(...)  ! OMESH writes
END SUBROUTINE
```

**Categorize each section**:
- OMESH reads → Phase 1 (preprocessing)
- Local cell-by-cell processing → Phase 2 (parallel kernel)
- OMESH writes / cross-mesh coupling → Phase 3 (finalization)

### Step 2: Convert Callees to Thread-Safe

All routines called in Phase 2 must be thread-safe. Use appropriate patterns:

#### Pattern A: Index-Based Access (Small/Medium Routines)

For routines with limited array access, use explicit indexing:

```fortran
! Before (module-level pointers)
FUNCTION HEAT_TRANSFER_COEFFICIENT(NM, ...)
  CALL POINT_TO_MESH(NM)
  WC => WALL(WALL_INDEX)
  B1 => BOUNDARY_PROP1(WC%B1_INDEX)
  B1%TMP_F = ...
END FUNCTION

! After (index-based, thread-safe)
FUNCTION HEAT_TRANSFER_COEFFICIENT(M, T, ...)
  TYPE(MESH_TYPE), INTENT(INOUT) :: M
  INTEGER :: B1_INDEX, B2_INDEX

  B1_INDEX = M%WALL(WALL_INDEX)%B1_INDEX
  M%BOUNDARY_PROP1(B1_INDEX)%TMP_F = ...  ! Direct array access
END FUNCTION
```

**Pros**: Simple, no language limitations
**Cons**: Verbose for large routines with many array accesses

#### Pattern B: Pointer-Based Access (Large Routines)

For routines with extensive array access, use local mesh pointer:

```fortran
! Key insight: MESHES is declared with TARGET attribute in MESH_VARIABLES
! TYPE(MESH_TYPE), SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: MESHES

SUBROUTINE SURFACE_HEAT_TRANSFER(NM, PREDICTOR_FLAG, T, ...)
  INTEGER, INTENT(IN) :: NM
  LOGICAL, INTENT(IN) :: PREDICTOR_FLAG
  TYPE(MESH_TYPE), POINTER :: M  ! ← POINTER, not INTENT(INOUT)
  REAL(EB), POINTER, DIMENSION(:,:,:) :: UU, VV, WW, RHOP

  M => MESHES(NM)  ! Works because MESHES has TARGET attribute

  ! Conditional pointer setup for predictor/corrector
  IF (PREDICTOR_FLAG) THEN
    UU => M%US; VV => M%VS; WW => M%WS; RHOP => M%RHOS
  ELSE
    UU => M%U; VV => M%V; WW => M%W; RHOP => M%RHO
  ENDIF

  ! Use pointers as before
  UN = UU(II,JJ,KK)
  VN = VV(II,JJ,KK)
END SUBROUTINE
```

**Critical**: Must use `TYPE(MESH_TYPE), POINTER :: M`, NOT `INTENT(INOUT)`. INTENT prevents pointer assignment to components.

**Pros**: Minimal changes to existing code, handles predictor/corrector elegantly
**Cons**: Requires MESHES to have TARGET attribute (FDS already has this)

### Step 3: Extract Preprocessing Routine

Create a new subroutine for Phase 1 operations:

```fortran
SUBROUTINE WALL_BC_PREPROCESSING(NM, T, DT_BC, CALL_HT_1D)
  INTEGER, INTENT(IN) :: NM
  REAL(EB), INTENT(IN) :: T, DT_BC
  LOGICAL, INTENT(IN) :: CALL_HT_1D

  TYPE(MESH_TYPE), POINTER :: M
  M => MESHES(NM)

  ! Cross-mesh reads (OMESH access allowed here)
  DO IW = 1, M%N_EXTERNAL_WALL_CELLS
    CALL ASSIGN_GHOST_VALUE(IW, ...)  ! Reads OMESH
  END DO

  ! Common setup for all cells
  DO IW = 1, M%N_EXTERNAL_WALL_CELLS + M%N_INTERNAL_WALL_CELLS
    CALL NEAR_SURFACE_GAS_VARIABLES_KERNEL(M, ...)
    IF (CALL_HT_1D .AND. SF%THERMAL_BC_INDEX == THERMALLY_THICK) THEN
      HEAT_COEF = HEAT_TRANSFER_COEFFICIENT(M, ...)
    ENDIF
  END DO
END SUBROUTINE
```

### Step 4: Extract Parallel Kernel

Create a new RECURSIVE kernel for Phase 2:

```fortran
RECURSIVE SUBROUTINE WALL_BC_PROCESS_CELLS_KERNEL(M, NM, PREDICTOR_FLAG, T, DT, DT_BC, CALL_HT_1D)
  TYPE(MESH_TYPE), POINTER :: M
  INTEGER, INTENT(IN) :: NM
  LOGICAL, INTENT(IN) :: PREDICTOR_FLAG, CALL_HT_1D
  REAL(EB), INTENT(IN) :: T, DT, DT_BC

  ! Process cells WITHOUT cross-mesh dependencies
  DO IW = 1, M%N_EXTERNAL_WALL_CELLS + M%N_INTERNAL_WALL_CELLS
    WC => M%WALL(IW)

    ! Skip cells requiring sequential processing
    IF (WC%HAS_INTERPOLATED_BC) CYCLE  ! Handled in preprocessing
    IF (WC%HAS_BACK_MESH) CYCLE        ! Handled in finalization

    ! Thread-safe processing (all callees are thread-safe)
    CALL SURFACE_HEAT_TRANSFER(NM, PREDICTOR_FLAG, ...)
    CALL CALCULATE_ZZ_F(NM, PREDICTOR_FLAG, ...)

    ! Additional cell types (CFACE, particles)
    ! ...
  END DO
END SUBROUTINE
```

**Key**: Mark as RECURSIVE for thread safety. Skip cells with cross-mesh flags.

### Step 5: Extract Finalization Routine

Create a new subroutine for Phase 3:

```fortran
SUBROUTINE WALL_BC_FINALIZE(NM, T, DT_BC, CALL_HT_1D)
  INTEGER, INTENT(IN) :: NM
  REAL(EB), INTENT(IN) :: T, DT_BC
  LOGICAL, INTENT(IN) :: CALL_HT_1D

  TYPE(MESH_TYPE), POINTER :: M
  M => MESHES(NM)

  ! Handle cross-mesh coupling (BACK_MESH)
  DO IW = 1, M%N_EXTERNAL_WALL_CELLS + M%N_INTERNAL_WALL_CELLS
    WC => M%WALL(IW)
    IF (WC%HAS_BACK_MESH) THEN
      CALL SOLID_HEAT_TRANSFER(..., BACK_MESH=...)  ! Cross-mesh access allowed
    ENDIF
  END DO

  ! Handle thin wall lateral heat transfer (all thin walls)
  DO IW = 1, M%N_THIN_WALL_CELLS
    CALL SOLID_HEAT_TRANSFER(...)
  END DO

  ! Particle off-gassing (OMESH writes)
  IF (CORRECTOR .AND. CALL_HT_1D) THEN
    CALL DEPOSIT_PARTICLE_MASS(NM, T, DT_BC)  ! Updates neighboring meshes
  ENDIF
END SUBROUTINE
```

### Step 6: Create C Wrappers

Add ISO_C_BINDING wrappers in `fds_c_interface.f90`:

```fortran
! Preprocessing wrapper
SUBROUTINE C_FDS_WALL_BC_PREPROCESSING(NM, T, DT_BC, CALL_HT_1D) &
    BIND(C, NAME="fds_wall_bc_preprocessing")
    USE WALL_ROUTINES, ONLY: WALL_BC_PREPROCESSING
    INTEGER(C_INT), VALUE :: NM, CALL_HT_1D
    REAL(C_DOUBLE), VALUE :: T, DT_BC
    LOGICAL :: HT_1D_FLAG
    HT_1D_FLAG = (CALL_HT_1D /= 0)
    CALL WALL_BC_PREPROCESSING(NM, T, DT_BC, HT_1D_FLAG)
END SUBROUTINE

! Kernel wrapper (RECURSIVE for thread safety)
RECURSIVE SUBROUTINE C_FDS_WALL_BC_PROCESS_CELLS_KERNEL(NM, T, DT, DT_BC, CALL_HT_1D) &
    BIND(C, NAME="fds_wall_bc_process_cells_kernel")
    USE WALL_ROUTINES, ONLY: WALL_BC_PROCESS_CELLS_KERNEL
    INTEGER(C_INT), VALUE :: NM, CALL_HT_1D
    REAL(C_DOUBLE), VALUE :: T, DT, DT_BC
    LOGICAL :: PREDICTOR_FLAG, HT_1D_FLAG
    PREDICTOR_FLAG = PREDICTOR  ! From global state
    HT_1D_FLAG = (CALL_HT_1D /= 0)
    CALL WALL_BC_PROCESS_CELLS_KERNEL(MESHES(NM), NM, PREDICTOR_FLAG, T, DT, DT_BC, HT_1D_FLAG)
END SUBROUTINE

! Finalization wrapper
SUBROUTINE C_FDS_WALL_BC_FINALIZE(NM, T, DT_BC, CALL_HT_1D) &
    BIND(C, NAME="fds_wall_bc_finalize")
    ! Similar structure
END SUBROUTINE
```

### Step 7: Create Helper Functions for Global State

For parameters computed from global state, create helper functions:

```fortran
! Compute DT_BC from global BC_CLOCK
FUNCTION C_FDS_COMPUTE_WALL_BC_DT_BC(T) RESULT(DT_BC_OUT) &
    BIND(C, NAME="fds_compute_wall_bc_dt_bc")
    USE MESH_VARIABLES, ONLY: BC_CLOCK
    REAL(C_DOUBLE), VALUE :: T
    REAL(C_DOUBLE) :: DT_BC_OUT
    DT_BC_OUT = T - BC_CLOCK
END FUNCTION

! Check if 1-D heat transfer should be called
FUNCTION C_FDS_CHECK_CALL_HT_1D() RESULT(CALL_HT_1D_OUT) &
    BIND(C, NAME="fds_check_call_ht_1d")
    USE MESH_VARIABLES, ONLY: WALL_COUNTER, WALL_INCREMENT
    INTEGER(C_INT) :: CALL_HT_1D_OUT
    IF (.NOT.INITIALIZATION_PHASE .AND. CORRECTOR) THEN
        IF (WALL_COUNTER == WALL_INCREMENT) THEN
            CALL_HT_1D_OUT = 1
        ELSE
            CALL_HT_1D_OUT = 0
        ENDIF
    ELSE
        CALL_HT_1D_OUT = 0
    ENDIF
END FUNCTION

! Update global BC_CLOCK
SUBROUTINE C_FDS_UPDATE_BC_CLOCK(T) BIND(C, NAME="fds_update_bc_clock")
    USE MESH_VARIABLES, ONLY: BC_CLOCK, HT_3D_SWEEP_DIRECTION
    REAL(C_DOUBLE), VALUE :: T
    BC_CLOCK = T
    HT_3D_SWEEP_DIRECTION = HT_3D_SWEEP_DIRECTION + 1
    IF (HT_3D_SWEEP_DIRECTION > 3) HT_3D_SWEEP_DIRECTION = 1
END SUBROUTINE
```

**Rationale**: Computing global parameters in the orchestrator (C++ side) avoids redundant computation and keeps state management clear.

### Step 8: Create Work Token Data Structure

Define a work token in `data/wallbc_data.h`:

```cpp
struct WallBCWork {
    int nm;              // Mesh index
    double t;            // Current time
    double dt;           // Time step
    double dt_bc;        // Boundary condition time step (computed globally)
    int call_ht_1d;      // 1-D heat transfer flag (computed globally)
    std::shared_ptr<MeshData> originalMeshData;  // For passing through

    WallBCWork(int nm_, double t_, double dt_, double dt_bc_, int call_ht_1d_,
               std::shared_ptr<MeshData> md)
        : nm(nm_), t(t_), dt(dt_), dt_bc(dt_bc_), call_ht_1d(call_ht_1d_),
          originalMeshData(md) {}
};
```

### Step 9: Create Hedgehog Orchestrator

Implement orchestrator in `state/wallbc_state.h`:

```cpp
class WallBCOrchestrator : public hh::AbstractState<1, MeshData, WallBCWork> {
public:
    explicit WallBCOrchestrator(int nmeshes)
        : hh::AbstractState<1, MeshData, WallBCWork>(),
          nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);

        if (static_cast<int>(collected_.size()) == nmeshes_) {
            // Compute global parameters once
            double dt_bc = fds_compute_wall_bc_dt_bc(collected_[0]->t);
            int call_ht_1d = fds_check_call_ht_1d();

            // Update global state if needed
            if (call_ht_1d) {
                fds_update_bc_clock(collected_[0]->t);
            }

            // Sequential preprocessing for each mesh
            for (auto &md : collected_) {
                fds_wall_bc_preprocessing(md->nm, md->t, dt_bc, call_ht_1d);
            }

            // Dispatch parallel work tokens
            for (auto &md : collected_) {
                auto work = std::make_shared<WallBCWork>(
                    md->nm, md->t, md->dt, dt_bc, call_ht_1d, md);
                this->addResult(work);
            }

            collected_.clear();
            collected_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};
```

### Step 10: Create Hedgehog Kernel Task

Implement kernel task in `task/wallbc_kernel_task.h`:

```cpp
class WallBCKernelTask : public hh::AbstractTask<1, WallBCWork, WallBCWork> {
public:
    explicit WallBCKernelTask(size_t numThreads)
        : hh::AbstractTask<1, WallBCWork, WallBCWork>("WallBCKernel", numThreads) {}

    void execute(std::shared_ptr<WallBCWork> work) override {
        fds_wall_bc_process_cells_kernel(
            work->nm, work->t, work->dt, work->dt_bc, work->call_ht_1d);
        this->addResult(work);
    }

    std::shared_ptr<hh::AbstractTask<1, WallBCWork, WallBCWork>> copy() override {
        return std::make_shared<WallBCKernelTask>(this->numberThreads());
    }
};
```

### Step 11: Create Hedgehog Collector

Implement collector in `state/wallbc_state.h`:

```cpp
class WallBCCollector : public hh::AbstractState<1, WallBCWork, MeshData> {
public:
    explicit WallBCCollector(int nmeshes)
        : hh::AbstractState<1, WallBCWork, MeshData>(),
          nmeshes_(nmeshes) {
        results_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<WallBCWork> work) override {
        results_.push_back(work);

        if (static_cast<int>(results_.size()) == nmeshes_) {
            // Sort by mesh index for deterministic ordering
            std::sort(results_.begin(), results_.end(),
                      [](const auto &a, const auto &b) {
                          return a->nm < b->nm;
                      });

            // Sequential finalization for each mesh
            for (auto &w : results_) {
                fds_wall_bc_finalize(w->nm, w->t, w->dt_bc, w->call_ht_1d);
            }

            // Emit MeshData tokens
            for (auto &w : results_) {
                this->addResult(w->originalMeshData);
            }

            results_.clear();
            results_.reserve(nmeshes_);
        }
    }

private:
    int nmeshes_;
    std::vector<std::shared_ptr<WallBCWork>> results_;
};
```

### Step 12: Wire into Main Graph

Add to `fds_graph.h`:

```cpp
// Create WallBC sub-graph components
auto wallBCOrchSM = std::make_shared<hh::StateManager<1, MeshData, WallBCWork>>(
    std::make_shared<WallBCOrchestrator>(nmeshes), "WallBCOrch");
auto wallBCKernelTask = std::make_shared<WallBCKernelTask>(kernelThreads);
auto wallBCCollectorSM = std::make_shared<hh::StateManager<1, WallBCWork, MeshData>>(
    std::make_shared<WallBCCollector>(nmeshes), "WallBCCollector");

// Wire sub-graph (replace sequential task)
graph->edges(prevNode, wallBCOrchSM);
graph->edges(wallBCOrchSM, wallBCKernelTask);
graph->edges(wallBCKernelTask, wallBCCollectorSM);
graph->edges(wallBCCollectorSM, nextNode);
```

### Step 13: Testing Strategy

1. **Build**: `cmake --build build_hh --target fds -j$(nproc)`

2. **1-mesh test** (sequential baseline):
   ```bash
   cd test_cases
   python3 run_tests.py -t dancing_eddies_1mesh -v
   ```

3. **Multi-mesh tests** (parallel execution):
   ```bash
   python3 run_tests.py -v  # All tests (1-5 meshes)
   ```

4. **Verification**:
   - All outputs must be byte-identical to gold files
   - Check _devc.csv and _hrr.csv for deterministic ordering
   - Verify no runtime errors or crashes

## Common Patterns and Tips

### Global Parameter Computation

**Pattern**: Compute once in orchestrator, pass to all work tokens

```cpp
// Orchestrator computes global parameters
double dt_bc = fds_compute_wall_bc_dt_bc(collected_[0]->t);

// All work tokens receive the same value
for (auto &md : collected_) {
    auto work = std::make_shared<WallBCWork>(..., dt_bc, ...);
}
```

**Avoid**: Computing global parameters in each kernel (redundant work, race conditions)

### Flag-Based Cell Filtering

**Pattern**: Use flags to separate parallel and sequential processing

```fortran
! In preprocessing: mark cells needing sequential processing
IF (BOUNDARY_TYPE == INTERPOLATED_BC) THEN
    WC%HAS_INTERPOLATED_BC = .TRUE.
ENDIF

! In kernel: skip flagged cells
IF (WC%HAS_INTERPOLATED_BC) CYCLE  ! Sequential only
IF (WC%HAS_BACK_MESH) CYCLE        ! Sequential only
```

### Sorting for Determinism

**Pattern**: Always sort results before emitting

```cpp
std::sort(results_.begin(), results_.end(),
          [](const auto &a, const auto &b) { return a->nm < b->nm; });
```

**Critical**: Hedgehog executes tasks in non-deterministic order. Sorting ensures reproducible output.

### Module-Level State Access

**Pattern**: Read global state in wrappers, pass as arguments

```fortran
! In C wrapper
PREDICTOR_FLAG = PREDICTOR  ! Read from global module
CALL KERNEL(..., PREDICTOR_FLAG)  ! Pass as argument

! In kernel
SUBROUTINE KERNEL(..., PREDICTOR_FLAG)
  LOGICAL, INTENT(IN) :: PREDICTOR_FLAG  ! No global access
```

**Avoid**: Reading module-level state inside parallel kernels (race conditions)

## Performance Expectations

For routines with 80-90% parallelizable work:

| Meshes | Threads | Expected Speedup |
|--------|---------|------------------|
| 1      | 1       | 1.0× (baseline)  |
| 4      | 4       | ~8-9× (super-linear, cache effects) |
| 8      | 8       | ~15-18×          |

**Sequential overhead**: 10-20% for preprocessing + finalization

## Troubleshooting

### Byte-Identical Verification Fails

**Cause**: Non-deterministic execution order
**Fix**: Ensure collector sorts results before finalization

### Runtime Crashes in Kernel

**Cause**: OMESH access in kernel, or non-thread-safe callee
**Fix**: Move OMESH access to preprocessing/finalization, verify callees are thread-safe

### Performance Below Expected

**Cause**: Too much work in preprocessing/finalization
**Fix**: Profile with Hedgehog dot file, move more work to parallel kernel

## Reference Implementation

**Complete example**: WallBC sub-graph
- Documentation: `docs/WALL_BC_PARALLELIZATION_PLAN.md`
- Test report: `test_cases/WALLBC_TEST_REPORT.md`
- Files modified: See commit `bf2163ee75`

## Summary Checklist

- [ ] Analyze routine: identify phases (OMESH reads → local → OMESH writes)
- [ ] Convert callees to thread-safe (index-based or pointer-based)
- [ ] Extract preprocessing routine (Phase 1, OMESH reads allowed)
- [ ] Extract parallel kernel (Phase 2, RECURSIVE, no OMESH)
- [ ] Extract finalization routine (Phase 3, OMESH writes allowed)
- [ ] Create C wrappers (preprocessing, kernel RECURSIVE, finalization)
- [ ] Create helper functions for global parameters
- [ ] Create work token data structure
- [ ] Implement orchestrator (collect → compute globals → preprocess → dispatch)
- [ ] Implement kernel task (call kernel, pass through work token)
- [ ] Implement collector (gather → sort → finalize → emit)
- [ ] Wire into main graph (replace sequential task)
- [ ] Test: 1-mesh byte-identical, multi-mesh byte-identical
- [ ] Document: method, test results, lessons learned
