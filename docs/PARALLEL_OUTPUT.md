# Parallel Output Plan — Final

## Goal

Eliminate the dump phase as a bottleneck by:
1. Reorganizing dump code for structural clarity
2. Removing redundant file open/close syscalls (persistent handles)
3. Making dump routines thread-safe (remove POINT_TO_MESH dependency)
4. Parallelizing per-mesh dump I/O across Hedgehog threads

## Current Architecture

### Hedgehog Integration (timestep_state.h)

`TimestepDumpState` is a monolithic state that runs after the corrector barrier:

```
BarrierData (all meshes collected)
  │
  ▼
TimestepDumpState:
  Phase 1: fds_set_diagnostics, fds_exchange_global_outputs, fds_update_controls  [global, sequential]
  Phase 2: for each mesh: fds_dump_mesh_outputs(t, dt, nm)                        [per-mesh, SEQUENTIAL ← bottleneck]
  Phase 3: fds_dump_global_outputs, fds_write_strings, fds_write_diagnostics       [global, sequential]
  Phase 4: termination decision, DT adjustment
```

### Fortran Dump Dispatch (dump.f90:88-218)

`DUMP_MESH_OUTPUTS(T, DT, NM)` mixes three concerns into one routine:
1. **Scheduling** — clock-based "should I dump now?" checks
2. **Context setup** — `POINT_TO_MESH(NM)` sets 400+ module-level pointer aliases
3. **Execution** — dispatches to 11 sub-routines that interleave computation and I/O

| Routine | Files/mesh/call | Open/Close per call? | Source |
|---------|----------------|---------------------|--------|
| `DUMP_PART` | 2 (.prt5 + .prt5.bnd) | YES | dump.f90:4134 |
| `DUMP_ISOF` | 2×N_ISOF (.iso + .viso) | YES | dump.f90:4256 |
| `DUMP_SMOKE3D` → `SMOKE3D_TO_FILE` | 3×N_SMOKE3D (.s3d + .s3d.sz + .s3dd) | YES | smvv.f90:924 |
| `DUMP_SLCF` (slices) | 3×N_SLCF (.sf + .sf.bnd + .sf.rle) | YES (11 OPEN, 10 CLOSE) | dump.f90:5968 |
| `DUMP_SLCF` (3D slices) | same pattern | YES | dump.f90:5968 |
| `DUMP_SLCF` (Plot3D) | 2 (.q + .q.bnd) | YES (STATUS='REPLACE') | dump.f90:6074 |
| `DUMP_BNDF` | 2×N_BNDF (.bf + .bf.bnd) + CC_IBM variants | YES | dump.f90:10582 |
| `DUMP_PROF` | 1 (.csv) | YES | dump.f90:10022 |
| `DUMP_UVW` | 1 (.csv, unique filename) | YES (STATUS='REPLACE') | dump.f90 |
| `DUMP_TMP` | 1 (.csv, unique filename) | YES (STATUS='REPLACE') | dump.f90 |
| `DUMP_SPEC` | 1 (.csv, unique filename) | YES (STATUS='REPLACE') | dump.f90 |

**Worst case per mesh per dump step**: ~20+ OPEN + 20+ CLOSE syscalls (if all output types fire).

### Internal Structure of Dump Routines

Each dump sub-routine currently interleaves computation and I/O:

```
DUMP_SLCF (543 lines, largest):
  OPEN file
  DO N = 1, N_SLCF              ← loop over slice definitions
    DO K/J/I = ...               ← loop over cells in slice
      val = GAS_PHASE_OUTPUT()   ← computation (pure, 0 I/O)
      QQ(I,J,K) = val            ← buffer
    END DO
    WRITE(LU) QQ                 ← file I/O
    OPEN/WRITE/CLOSE bounds      ← file I/O
  END DO
  CLOSE file
```

However, the hub routines are already clean:
- **GAS_PHASE_OUTPUT** (1432 lines): pure computation, 0 I/O, takes NM parameter
- **SOLID_PHASE_OUTPUT** (826 lines): pure computation, 0 I/O, takes NM parameter
- **UPDATE_GLOBAL_OUTPUTS** (20 lines): pure computation, 0 I/O, calls POINT_TO_MESH

### Global Files (already optimized)

Global files (LU_HRR, LU_MASS, LU_DEVC, LU_CTRL, LU_HVAC, LU_STEPS) are
opened once during `INITIALIZE_GLOBAL_DUMPS` and kept open with periodic
`FLUSH_GLOBAL_BUFFERS`. No changes needed.

### Thread-Safety Blocker

Only two call sites invoke `POINT_TO_MESH(NM)`:
1. `DUMP_MESH_OUTPUTS` (dump.f90:99)
2. `UPDATE_GLOBAL_OUTPUTS` (dump.f90:74)

The sub-routines themselves do NOT call POINT_TO_MESH — they rely on the
module aliases already being set by the caller. They access data through
module-level pointers like `U`, `V`, `W`, `RHO`, `CELL`, `WALL`, `WORK1`,
`QQ`, `IBP1`, `XC`, etc.

---

## Implementation Plan

### Phase 1: Reorganization

**Goal**: Restructure dump code so that scheduling, file management,
computation, and I/O are separated concerns. This makes every subsequent
phase cleaner and less error-prone.

**Scope**: `dump.f90`, `smvv.f90`

#### 1a. Decouple scheduling from execution in DUMP_MESH_OUTPUTS

Currently `DUMP_MESH_OUTPUTS` is a 130-line routine that mixes clock checks,
counter advancement, and routine dispatch. Split into two layers:

```fortran
! Layer 1: Pure scheduling — what needs dumping this timestep?
! Returns a bitmask or set of flags.
SUBROUTINE CHECK_DUMP_SCHEDULE(T, NM, DO_PART, DO_ISOF, DO_SM3D, DO_SLCF, &
                                DO_SL3D, DO_BNDF, DO_PL3D, DO_PROF)
  LOGICAL, INTENT(OUT) :: DO_PART, DO_ISOF, DO_SM3D, DO_SLCF, &
                           DO_SL3D, DO_BNDF, DO_PL3D, DO_PROF

  DO_PART = (T >= PART_CLOCK(PART_COUNTER(NM)) .AND. PARTICLE_FILE)
  DO_ISOF = (T >= ISOF_CLOCK(ISOF_COUNTER(NM)))
  DO_SM3D = (T >= SM3D_CLOCK(SM3D_COUNTER(NM)) .AND. SMOKE3D)
  ! ... etc ...
END SUBROUTINE

! Layer 2: Advance counters (called after execution)
SUBROUTINE ADVANCE_DUMP_COUNTERS(T, NM, DO_PART, DO_ISOF, ...)
  IF (DO_PART) THEN
    DO WHILE(PART_COUNTER(NM) < SIZE(PART_CLOCK)-1)
      PART_COUNTER(NM) = PART_COUNTER(NM) + 1
      IF (PART_CLOCK(PART_COUNTER(NM)) >= T) EXIT
    END DO
  END IF
  ! ... etc ...
END SUBROUTINE
```

**Why**: The scheduling decision is needed by the Hedgehog graph to know
whether to dispatch dump work at all (skip entirely on non-dump timesteps).
The counter advancement must happen after execution. Mixing them in one
routine prevents this separation.

#### 1b. Separate computation from I/O in dump routines

Within each dump routine, split the compute-and-buffer step from the
write-to-file step. The compute step fills output buffers; the write step
flushes them to disk. Both stay in the same subroutine for now but are
clearly delineated sections.

**DUMP_SLCF** (the largest, 543 lines):

Currently interleaves per-slice: compute cells → write → compute bounds → write bounds.
Restructure as:

```fortran
SUBROUTINE DUMP_SLCF(T, DT, NM, IFRMT)
  ! --- COMPUTE PHASE ---
  ! Build solid mask B, normalization S (uses WORK arrays)
  CALL PREPARE_SLCF_MASK(NM)

  SLICE_LOOP: DO N = 1, N_SLCF
    ! Fill output buffer QQ via GAS_PHASE_OUTPUT
    CALL COMPUTE_SLCF_QUANTITY(T, DT, NM, N, QQ, SLICE_MIN, SLICE_MAX)
  END DO SLICE_LOOP

  ! --- WRITE PHASE ---
  SLICE_LOOP_WRITE: DO N = 1, N_SLCF
    CALL WRITE_SLCF_DATA(NM, N, QQ, SLICE_MIN, SLICE_MAX)
  END DO SLICE_LOOP_WRITE
END SUBROUTINE
```

In practice, the two loops can remain fused (compute+write per slice) since
each slice uses the same WORK/QQ arrays and we'd need per-slice buffers to
fully separate them. The key change is extracting the compute and write
logic into clearly named internal subroutines (via CONTAINS), so each can
later be called independently.

**DUMP_BNDF**, **DUMP_ISOF**, **DUMP_PART**, **DUMP_SMOKE3D**: Same pattern.
Extract `COMPUTE_*` and `WRITE_*` internal procedures.

**Why**: Once compute and write are separate:
- Phase 2 (thread safety) only needs to convert the compute functions
- Phase 3 (parallelism) can run compute in parallel and serialize writes
  if needed (though for per-mesh files, writes can also be parallel)
- The code becomes easier to read: "what data is computed" vs "where it goes"

#### 1c. Group file lifecycle management

Currently file initialization (header writing), runtime opens, and implicit
close-at-exit are scattered across three places. Consolidate:

**Per output type, create paired routines**:

```fortran
! Called from INITIALIZE_MESH_DUMPS
SUBROUTINE OPEN_SLCF_FILES(NM, RESTART)
  ! Fresh run: create with STATUS='REPLACE', write headers
  ! Restart: open with STATUS='OLD', POSITION='APPEND'
END SUBROUTINE

! Called from finalization
SUBROUTINE CLOSE_SLCF_FILES(NM)
  DO N = 1, MESHES(NM)%N_SLCF
    CLOSE(LU_SLCF(N,NM))
    CLOSE(LU_SLCF(N+N_SLCF_MAX,NM))
    CLOSE(LU_SLCF(N+2*N_SLCF_MAX,NM))
  END DO
END SUBROUTINE
```

Do this for SLCF, BNDF, ISOF, PART, SMOKE3D, PL3D.

Then the master routines become:
```fortran
SUBROUTINE OPEN_ALL_MESH_OUTPUT_FILES(NM, RESTART)
  CALL OPEN_SLCF_FILES(NM, RESTART)
  CALL OPEN_BNDF_FILES(NM, RESTART)
  CALL OPEN_ISOF_FILES(NM, RESTART)
  CALL OPEN_PART_FILES(NM, RESTART)
  CALL OPEN_SMOKE3D_FILES(NM, RESTART)
END SUBROUTINE

SUBROUTINE CLOSE_ALL_MESH_OUTPUT_FILES(NM)
  CALL CLOSE_SLCF_FILES(NM)
  CALL CLOSE_BNDF_FILES(NM)
  ! ... etc ...
END SUBROUTINE
```

**Why**: This makes Phase 2 trivial — persistent handles just means removing
the OPEN/CLOSE from dump routines and relying on the lifecycle routines.
It also makes the file state explicit and auditable.

#### 1d. Isolate GAS_PHASE_OUTPUT / SOLID_PHASE_OUTPUT alias dependencies

These hub functions are already pure computation, but they access mesh data
through module-level aliases set by POINT_TO_MESH. Document exactly which
aliases each function uses, as preparation for Phase 3's thread-safe conversion.

Create a comment block at the top of each function:

```fortran
! Module aliases used (set by POINT_TO_MESH):
!   From MESH_VARIABLES: U, V, W, US, VS, WS, RHO, TMP, ZZ, Q, QR,
!     UII, D, DS, H, HS, KRES, MU, PBAR, D_PBAR_DT, KAPPA_GAS,
!     FVX, FVY, FVZ, CELL, CELL_INDEX, CELL_COUNT, EDGE,
!     IBP1, JBP1, KBP1, IBAR, JBAR, KBAR, XC, YC, ZC, X, Y, Z,
!     DX, DY, DZ, DXN, DYN, DZN, RDXN, RDYN, RDZN
!   From CC_SCALARS (if CC_IBM): CUT_FACE, FCVAR, CCVAR, CC_VGSC
```

Additionally, note that GAS_PHASE_OUTPUT already takes `NM` as a parameter
and accesses `MESHES(NM)` in a few places (e.g., for CC_IBM paths). This
means it can be converted to use `M => MESHES(NM)` throughout in Phase 3.

#### 1e. Testing

- All changes are pure refactoring — zero behavioral change
- Run all 18 custom tests + verification suite
- Verify identical output files (binary diff)

**Risk**: Low. Pure structural refactoring, no logic changes.

---

### Phase 2: Persistent File Handles

**Goal**: Eliminate open/close syscalls. Open files once, keep them open,
close at finalization.

**Scope**: `dump.f90`, `smvv.f90`

**Prerequisite**: Phase 1 (file lifecycle routines already grouped)

#### 2a. Remove the negative-LU sign convention for ISOF

Currently (dump.f90:394-397):
```fortran
LU_ISOF(N,NM) = -GET_FILE_NUMBER()
IF (RESTART) LU_ISOF(N,NM) = ABS(LU_ISOF(N,NM))
```
The negative sign tells `DUMP_ISOF` whether this is the first call (create
file with REPLACE) or subsequent (append). With persistent handles, files
are always already open, so this convention is unnecessary.

**Action**: Use positive LU numbers directly. The OPEN_ISOF_FILES routine
from Phase 1 handles both fresh and restart cases.

#### 2b. Keep files open after initialization

In the `OPEN_*_FILES` routines created in Phase 1, simply leave files open
after writing headers (remove the CLOSE calls that currently follow header
writes).

For restart runs: open with `STATUS='OLD', POSITION='APPEND'` and leave open.

This is now a trivial change because Phase 1 already grouped all opens/closes.

#### 2c. Remove OPEN/CLOSE from dump sub-routines

Strip all `OPEN(... STATUS='OLD', POSITION='APPEND')` and matching `CLOSE`
calls from the WRITE_* internal procedures created in Phase 1:
- `DUMP_PART` / `WRITE_PART_DATA`
- `DUMP_ISOF` / `WRITE_ISOF_DATA`
- `DUMP_SLCF` / `WRITE_SLCF_DATA`
- `DUMP_BNDF` / `WRITE_BNDF_DATA`
- `SMOKE3D_TO_FILE`

The write statements continue to use the same LU unit numbers — they're just
already open.

**Exception**: Plot3D uses `STATUS='REPLACE'` (overwrites each time). Use
`REWIND` + `WRITE` on persistent handles instead.

**Exception**: `DUMP_UVW`, `DUMP_TMP`, `DUMP_SPEC` create unique filenames
per dump step (e.g. `_uvw_t3_m1.csv`). These must remain open/write/close
since each dump produces a new file.

#### 2d. Bounds files special handling

Boundary bounds files (`.bf.bnd`, `.sf.bnd`) use read-modify-write:
```fortran
OPEN(...STATUS='OLD')
READ(LU) old_min, old_max
new_min = MIN(old_min, current_min)
CLOSE(LU)
OPEN(...STATUS='REPLACE')
WRITE(LU) new_min, new_max
CLOSE(LU)
```
With persistent handles, use `REWIND` before read/write cycles instead of
close+reopen with REPLACE.

#### 2e. Flush strategy

With persistent handles, data stays in OS buffers longer. Add FLUSH calls
in `FLUSH_GLOBAL_BUFFERS` (already called periodically) for per-mesh files:

```fortran
! Add to existing FLUSH_GLOBAL_BUFFERS or create FLUSH_MESH_BUFFERS
DO NM = LOWER_MESH_INDEX, UPPER_MESH_INDEX
  DO N = 1, MESHES(NM)%N_SLCF
    FLUSH(LU_SLCF(N,NM))
  END DO
  ! ... etc for BNDF, ISOF, PART, SMOKE3D ...
END DO
```

This ensures data reaches disk periodically for crash recovery.

#### 2f. Testing

- Run all 18 custom tests + verification suite
- Binary-compare output files (`.sf`, `.bf`, `.prt5`, `.iso`) against Phase 1 baseline
- Restart test: run 5 steps, restart, compare output continuity
- On NFS: verify no data loss from buffering

**Risk**: Low. File contents unchanged, only syscall pattern changes.

---

### Phase 3: Thread-Safe Dump Routines

**Goal**: Remove POINT_TO_MESH dependency so dump routines can run in parallel.

**Scope**: `dump.f90`, `smvv.f90`

**Prerequisite**: Phase 1 (compute/write already separated, alias audit done)

**Pattern**: Same approach used for all other Hedgehog kernels
(see MEMORY.md "Thread-Safe Fortran Patterns"):

```fortran
! Before (thread-unsafe — relies on module aliases set by POINT_TO_MESH):
SUBROUTINE DUMP_SLCF(T, DT, NM, IFRMT)
  ... uses U, V, W, RHO, CELL, WORK1, QQ, IBP1, XC etc. ...

! After (thread-safe — local pointer to mesh):
SUBROUTINE DUMP_SLCF(T, DT, NM, IFRMT)
  TYPE(MESH_TYPE), POINTER :: M
  M => MESHES(NM)
  ... uses M%U, M%V, M%W, M%RHO, M%CELL, M%WORK1, M%QQ, M%IBP1, M%XC etc. ...
```

#### 3a. Convert GAS_PHASE_OUTPUT and SOLID_PHASE_OUTPUT

These are the hub functions called by all dump routines. Convert first since
all dump routines depend on them.

**GAS_PHASE_OUTPUT** (1432 lines, 259 CASE branches):
- Already takes NM as a parameter
- Add `TYPE(MESH_TYPE), POINTER :: M` local variable, set `M => MESHES(NM)`
- Replace all module aliases with M% references
- The alias audit from Phase 1d provides the complete list

**SOLID_PHASE_OUTPUT** (826 lines, 60+ CASE branches):
- Already takes NM as a parameter
- Same M% conversion pattern

**Risk note**: GAS_PHASE_OUTPUT is also called by UPDATE_DEVICES_1 (device
output). After conversion, UPDATE_DEVICES_1 no longer needs POINT_TO_MESH
either — it already passes NM to GAS_PHASE_OUTPUT. Verify device output
is unaffected.

#### 3b. Convert dump sub-routines

Thanks to Phase 1b, each dump routine's computation is already in a
COMPUTE_* internal procedure. Convert these to use M%:

Priority order (by complexity):

1. **DUMP_SLCF** (543 lines) — uses WORK1→B, WORK2→S, WORK3→QUANTITY, QQ,
   CELL, CELL_INDEX, CELL_COUNT, IBP1/JBP1/KBP1, XC/YC/ZC, X/Y/Z,
   DX/DY/DZ, IBAR/JBAR/KBAR. Calls GAS_PHASE_OUTPUT (already converted).
   Also accesses M2 => MESHES(NOM) for cross-mesh lookups — this is fine,
   each mesh's M pointer is local.

2. **DUMP_BNDF** (220 lines) — uses WALL, BOUNDARY_PROP1, BOUNDARY_COORD,
   BOUNDARY_ONE_D, CELL. Calls SOLID_PHASE_OUTPUT (already converted).

3. **DUMP_ISOF** (181 lines) — uses CELL, CELL_INDEX, GAS_PHASE_OUTPUT,
   WORK arrays. Already uses MESHES(NM) explicitly in some places.

4. **DUMP_PART** (114 lines) — uses LAGRANGIAN_PARTICLE, NLP, OMESH,
   PARTICLE_TAG. Note: LAGRANGIAN_PARTICLE_CLASS is a global array,
   not per-mesh — no M% needed for it.

5. **DUMP_SMOKE3D** / **SMOKE3D_TO_FILE** (120 lines) — uses GAS_PHASE_OUTPUT,
   CELL, WORK3, QQ.

6. **DUMP_PROF**, **DUMP_UVW**, **DUMP_TMP**, **DUMP_SPEC** — small routines.

#### 3c. Create DUMP_MESH_OUTPUTS_TS

New thread-safe entry point that replaces DUMP_MESH_OUTPUTS:

```fortran
SUBROUTINE DUMP_MESH_OUTPUTS_TS(T, DT, NM)
  ! No POINT_TO_MESH — each sub-routine uses M => MESHES(NM) internally
  ! Scheduling + dispatch (uses CHECK_DUMP_SCHEDULE from Phase 1)
  CALL CHECK_DUMP_SCHEDULE(T, NM, DO_PART, DO_ISOF, ...)

  IF (DO_PART) CALL DUMP_PART(T, NM)
  IF (DO_ISOF) CALL DUMP_ISOF(T, DT, NM)
  IF (DO_SM3D) CALL DUMP_SMOKE3D(T, DT, NM)
  IF (DO_SLCF) CALL DUMP_SLCF(T, DT, NM, 0)
  ! ... etc ...

  CALL ADVANCE_DUMP_COUNTERS(T, NM, DO_PART, DO_ISOF, ...)
END SUBROUTINE
```

#### 3d. Testing

- Run verification suite, compare against Phase 2 baseline
- Temporarily run with 1 thread to verify correctness before Phase 4
- Thread-safety verified by multi-thread execution in Phase 4

**Risk**: Medium. Large diff touching many routines. GAS_PHASE_OUTPUT is
the riskiest — 1432 lines of module alias replacements. Phase 1d's audit
mitigates this by documenting all aliases upfront.

---

### Phase 4: Parallel Dump via Hedgehog

**Goal**: Run per-mesh dump I/O in parallel using Hedgehog task threads.

**Scope**: `timestep_state.h`, `timestep_tasks.h`, `main_hh.cpp`

**Prerequisite**: Phase 3 (thread-safe dump routines)

#### 4a. Split TimestepDumpState into 3 parts

Replace the monolithic `TimestepDumpState` with a pipeline:

```
BarrierData
  │
  ▼
TimestepPreDumpState (state, 1 thread)
  │  fds_set_diagnostics
  │  fds_exchange_global_outputs
  │  fds_update_controls
  │  Emits individual MeshData tokens
  │
  ▼ (multicast to N mesh threads)
DumpMeshOutputsTask (task, N threads)
  │  fds_dump_mesh_outputs_ts(t, dt, nm)  ← thread-safe
  │  Each mesh writes to its own files
  │
  ▼ (barrier collect)
TimestepPostDumpState (state, 1 thread)
     fds_dump_global_outputs
     fds_write_strings
     fds_write_diagnostics
     fds_stop_check
     Termination decision + DT adjustment
     Emits BarrierData → TimestepLoopState
```

#### 4b. Thread count for DumpMeshOutputsTask

Use `nmeshes` threads (same as kernel tasks). Each thread handles one mesh
at a time. Dump I/O is I/O-bound, not CPU-bound, so having threads >= meshes
ensures no mesh waits.

On systems with slow I/O (NFS), having many threads writing simultaneously
could cause contention. Consider a configurable cap (e.g. `--dump-threads N`).

#### 4c. UPDATE_GLOBAL_OUTPUTS placement

Currently, `UPDATE_GLOBAL_OUTPUTS` runs per-mesh (calls POINT_TO_MESH)
and accumulates global quantities (HRR, mass) via UPDATE_HRR, UPDATE_MASS,
UPDATE_DEVICES_1. These routines update shared global arrays.

**Strategy**: Keep `UPDATE_GLOBAL_OUTPUTS` sequential in
`TimestepPreDumpState` (before parallel dump). It reads mesh data and
accumulates into global arrays — this must complete before DUMP_GLOBAL_OUTPUTS
writes the results. Since UPDATE_GLOBAL_OUTPUTS is pure computation (0 I/O),
it's fast and not worth parallelizing independently.

After Phase 3a converts GAS_PHASE_OUTPUT / SOLID_PHASE_OUTPUT to use M%
internally, UPDATE_GLOBAL_OUTPUTS no longer needs POINT_TO_MESH either.
However, UPDATE_HRR, UPDATE_MASS, UPDATE_DEVICES_1 themselves may still
use module aliases — these need the same M% treatment if we want to
parallelize them later (out of scope for this plan).

#### 4d. Graph wiring

```cpp
// Pre-dump state (barrier → mesh tokens)
auto preDumpState = std::make_shared<TimestepPreDumpState>(...);
auto preDumpSM = std::make_shared<TimestepPreDumpStateManager>(preDumpState, "PreDump");

// Parallel dump task
auto dumpTask = std::make_shared<DumpMeshOutputsTask>(nmeshes);  // N threads

// Post-dump barrier + finalization
auto postDumpCollector = makeBarrierSM(nmeshes, "PostDumpCollect", ...);
auto postDumpState = std::make_shared<TimestepPostDumpState>(...);

// Wiring
graph->edge<BarrierData>(correctorBarrier, preDumpSM);
graph->edge<MeshData>(preDumpSM, dumpTask);
graph->edge<MeshData>(dumpTask, postDumpCollector);
graph->edge<BarrierData>(postDumpCollector, postDumpState);
graph->edge<BarrierData>(postDumpState, timestepLoopSM);  // cycle back
```

#### 4e. Skip-dump optimization

On non-dump timesteps (no clocks fire), the entire dump pipeline is wasted
overhead. Use the scheduling separation from Phase 1a:

```cpp
void TimestepPreDumpState::execute(std::shared_ptr<BarrierData> data) {
    // ... global ops ...

    // Check if ANY mesh needs dumping this timestep
    bool anyDump = false;
    for (auto &md : data->meshes) {
        if (fds_check_any_dump(md->t, md->nm)) { anyDump = true; break; }
    }

    if (anyDump) {
        // Emit MeshData for parallel dump
        for (auto &md : data->meshes) this->addResult(md);
    } else {
        // Skip dump entirely — emit BarrierData directly to post-dump
        this->addResult(data);  // needs dual output type
    }
}
```

This avoids the barrier-collect-emit overhead on the ~90% of timesteps
where no output is needed.

#### 4f. Testing

- Verify bit-identical output with sequential dump (Phase 3 baseline)
- Profile wall-clock time for dump phase: compare 1 thread vs N threads
- Stress test with many meshes (16+) and frequent dump intervals
- Test skip-dump path with infrequent output intervals

**Risk**: Low (once Phase 3 is complete). The parallelism is embarrassingly
parallel — each mesh writes to independent files.

---

## Implementation Order & Dependencies

```
Phase 1: Reorganization
  ├── 1a: Split DUMP_MESH_OUTPUTS into schedule/dispatch/advance
  ├── 1b: Extract COMPUTE_*/WRITE_* within dump routines
  ├── 1c: Group OPEN_*/CLOSE_* file lifecycle routines
  ├── 1d: Audit and document alias dependencies
  └── 1e: Test (identical output)
           │
           ▼
Phase 2: Persistent File Handles
  ├── 2a: Remove negative-LU convention
  ├── 2b: Keep files open after initialization
  ├── 2c: Remove OPEN/CLOSE from WRITE_* routines
  ├── 2d: Bounds files REWIND handling
  ├── 2e: Add FLUSH strategy
  └── 2f: Test (identical output)
           │
           ▼
Phase 3: Thread-Safe Dump Routines
  ├── 3a: Convert GAS_PHASE_OUTPUT, SOLID_PHASE_OUTPUT to M%
  ├── 3b: Convert dump sub-routines (SLCF, BNDF, ISOF, PART, SMOKE3D)
  ├── 3c: Create DUMP_MESH_OUTPUTS_TS entry point
  └── 3d: Test (identical output, single-threaded)
           │
           ▼
Phase 4: Parallel Dump via Hedgehog
  ├── 4a: Split TimestepDumpState into Pre/Post
  ├── 4b: Configure dump thread count
  ├── 4c: Handle UPDATE_GLOBAL_OUTPUTS placement
  ├── 4d: Graph wiring
  ├── 4e: Skip-dump optimization
  └── 4f: Test & profile
```

Each phase produces identical output files. Each phase can be tested,
committed, and validated independently before proceeding.

## Files to Modify

| Phase | File | Changes |
|-------|------|---------|
| 1 | `Source/dump.f90` | Split DUMP_MESH_OUTPUTS; extract COMPUTE/WRITE internals; group OPEN/CLOSE per type; alias audit comments |
| 1 | `Source/smvv.f90` | Extract COMPUTE/WRITE in SMOKE3D_TO_FILE |
| 2 | `Source/dump.f90` | Remove OPEN/CLOSE from WRITE_* routines; keep OPEN_* files open; add CLOSE_ALL; add FLUSH |
| 2 | `Source/smvv.f90` | Remove OPEN/CLOSE from SMOKE3D_TO_FILE |
| 2 | `Source/hedgehog/fds_c_interface.f90` | Add C binding for CLOSE_ALL_MESH_OUTPUT_FILES |
| 2 | `Source/hedgehog/main_hh.cpp` | Call close_all_mesh_output_files after graph |
| 3 | `Source/dump.f90` | Convert GAS_PHASE_OUTPUT, SOLID_PHASE_OUTPUT + 6 dump routines to M% |
| 3 | `Source/smvv.f90` | Convert SMOKE3D_TO_FILE to M% |
| 3 | `Source/hedgehog/fds_c_interface.f90` | Add C binding for DUMP_MESH_OUTPUTS_TS |
| 4 | `Source/hedgehog/state/timestep_state.h` | Split TimestepDumpState into Pre/Post |
| 4 | `Source/hedgehog/task/timestep_tasks.h` | Update DumpMeshOutputsTask to use TS version |
| 4 | `Source/hedgehog/main_hh.cpp` | Rewire graph with parallel dump pipeline |

## Expected Performance Impact

**Phase 1** (reorganization):
- Zero performance change — pure refactoring
- Enables all subsequent optimizations

**Phase 2** (persistent handles):
- Eliminates ~20 OPEN + 20 CLOSE syscalls per mesh per dump step
- For 100 meshes, 1000 dump steps: eliminates ~4M syscalls
- Most impactful on NFS/network filesystems where open() latency is high
- Estimated improvement: 10-30% reduction in dump wall time

**Phase 4** (parallel dump):
- Dump I/O scales with number of threads (up to number of meshes)
- For N meshes: dump time ≈ max(single mesh dump time) instead of sum
- Estimated improvement: near-linear speedup up to I/O bandwidth saturation
- Combined with Phase 2: dump phase becomes negligible for most simulations

## Open Questions

1. **GEOM output** (`DUMP_GEOM`, dump.f90:10808): Uses unstructured geometry
   files. Needs investigation for persistent handles and thread safety.
   Lower priority since it's only relevant for CC_IBM cases.

2. **Restart files**: `DUMP_RESTART` uses STATUS='REPLACE' and is infrequent.
   Leave as-is (no persistent handle needed).

3. **File descriptor limits**: With many meshes and output types, persistent
   handles could hit OS `ulimit -n` limits. For 100 meshes with 20 files each
   = 2000 FDs. Default Linux limit is 1024. May need `ulimit -n 4096` or
   equivalent. Document this requirement.

4. **Flush strategy**: With persistent handles, data stays in OS buffers
   longer. Periodic FLUSH calls ensure data reaches disk for crash recovery.
   Frequency TBD — every dump step is safest but adds syscalls back;
   every N dump steps is a tradeoff.

5. **GAS_PHASE_OUTPUT cross-mesh access**: Some CASE branches in
   GAS_PHASE_OUTPUT access MESHES(NOM) for neighboring mesh data. This is
   safe for parallel execution (read-only access to other meshes), but needs
   verification that no branch writes to another mesh's arrays.
