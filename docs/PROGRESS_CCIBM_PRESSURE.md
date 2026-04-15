# CC_IBM Pressure Subgraph — Progress

## Goal

Enable the parallel pressure iteration subgraph (`buildPressureIterationSubgraph`)
for CC_IBM cases.  The CC_IBM gate in `fds_use_pressure_subgraph()` has been
removed, and CC_IBM calls (items 2-4) are integrated into pressure tasks.
All items complete — pressure subgraph fully supports CC_IBM (same-rank).

## Current State

The non-CC_IBM pressure pipeline is:

```
BaroclinicKernel → PreSolveExchange → PressureSolve → PostSolveExchange
    → VelocityError → ConvergenceCheck → (cycle or exit)
```

CC_IBM adds extra operations at 4 points in this pipeline.  All must be
integrated before the subgraph can be enabled.

## Work Items

### 1. Pre-loop exchange + GET_LINKED_FV

**Status**: ✅ DONE
**Complexity**: Low (~10 lines wiring)

```fortran
! main.f90:1423-1429 — runs ONCE before PRESSURE_ITERATION_LOOP
IF (CC_IBM) THEN
   CALL MESH_EXCHANGE(5)
   DO NM=...
      CALL GET_LINKED_FV(NM, DO_BAROCLINIC=.FALSE.)
   ENDDO
ENDIF
```

Links cut-face velocity fluxes across mesh boundaries before the iteration
loop starts.

**Plan**: Add to the upstream barrier where `fds_pressure_iteration_init()` +
`fds_pressure_iteration_increment()` are already called (in
`predJoinDivExchangeSM` / `predDivExchangeSM` for predictor,
`corrDivExchangeSM` for corrector).  This is already a barrier context, so
just append the call.

**Thread safety**: N/A — runs inside a barrier. GET_LINKED_FV is already thread-safe
(uses M% pattern, no POINT_TO_MESH).

**Implementation**: Added `fds_get_linked_fv(md->nm, 0)` per-mesh loop inside the
`if (useParallelPressure)` block of both CC_IBM barriers (predictor "WallDiv+DivExch"
and corrector "CorrDivExchange"), after `fds_pressure_iteration_init/increment`.
No pre-loop MESH_EXCHANGE(5) needed — GET_LINKED_FV reads only current mesh data;
the first pressure iteration's pre-solve exchange handles neighbor data.

**Files**:
- `fds_c_interface.f90`: C wrapper `fds_get_linked_fv` (line ~1294)
- `fds_fortran_interface.h`: declare wrapper
- `predictor_subgraph.h`: add call to CC_IBM barrier lambdas
- `corrector_subgraph.h`: same

### 2. CC_NO_FLUX in baroclinic phase

**Status**: ✅ DONE
**Complexity**: Medium (~150 lines Fortran)

Integrated into `BaroclinicKernelTask::execute()` (baroclinic_kernel_task.h:58-60).
Calls `fds_cc_no_flux(md->dt, md->nm, 1)` with FORCE_FLG=TRUE after baroclinic correction.

**Files**:
- `fds_c_interface.f90`: C wrapper `fds_cc_no_flux` (line 1273-1278)
- `fds_fortran_interface.h`: declaration (line 179)
- `task/baroclinic_kernel_task.h`: call site

### 3. CC_NO_FLUX in solve phase

**Status**: ✅ DONE (reuses item 2)

Integrated into `PressureSolveKernelTask::execute()` (pressure_iteration_tasks.h:63-65).
Calls `fds_cc_no_flux(md->dt, md->nm, 0)` with FORCE_FLG=FALSE after no_flux_kernel.

**Files**:
- `task/pressure_iteration_tasks.h`: call site

### 4. CC_COMPUTE_VELOCITY_ERROR

**Status**: ✅ DONE

Integrated into `VelocityErrorTask::execute()` (velocity_error_task.h:46-48).
Calls `fds_cc_compute_velocity_error(data->dt, data->nm)` after velocity error kernel.

**Files**:
- `fds_c_interface.f90`: C wrapper `fds_cc_compute_velocity_error` (line 1286-1292)
- `fds_fortran_interface.h`: declaration (line 181)
- `task/velocity_error_task.h`: call site

### 5. MESH_CC_EXCHANGE(5) — cut-face data exchange

**Status**: ✅ DONE (same-rank only; MPI cross-rank handled by barrier mode)

For same-rank (single MPI process), the CODE=5 CC exchange is simpler than
originally analyzed.  The same-rank path only copies `FN → FN_OMESH` on the
SENDER mesh's own cut-faces (H/HS read directly from neighbor mesh).

**Chosen approach — Option C: Pre-exchange kernel prep**:

Instead of modifying FluxExchangeTask or adding a barrier, set `FN_OMESH = FN`
in `BaroclinicKernelTask` BEFORE data enters the exchange cycle.  This works
because MeshDepsManager ensures all neighbors have arrived (and thus set their
FN_OMESH) before FluxExchangeTask runs for any mesh.

**Implementation**:
- `C_FDS_CC_EXCHANGE_PREPARE_FN(NM)` in `fds_c_interface.f90`: iterates over
  all neighbors via `N_NEIGHBORING_MESHES`, copies `CF%FN → CF%FN_OMESH` for
  boundary cut-faces.  Thread-safe (writes only to own mesh's CUT_FACE).
- Called in `BaroclinicKernelTask::doWork()` after `fds_cc_no_flux()`.
- MPI mode: redundant (barrier's `fds_mesh_exchange(5)` already handles CC),
  but harmless (negligible cost).

**Files**:
- `fds_c_interface.f90`: `C_FDS_CC_EXCHANGE_PREPARE_FN` routine
- `fds_fortran_interface.h`: `fds_cc_exchange_prepare_fn` declaration
- `task/baroclinic_kernel_task.h`: call site after `fds_cc_no_flux`

### 6. Enable the subgraph

**Status**: ✅ DONE (all items complete)
**Complexity**: Low (~5 lines)

The `IF (CC_IBM) RETURN` guard has been removed from `C_FDS_USE_PRESSURE_SUBGRAPH()`
(fds_c_interface.f90:145-158). Function now only checks: TUNNEL_PRECONDITIONER,
PRES_FLAG (FFT or ULMAT), and mesh count (≥2).

All prerequisite items (1-5) are complete.  Tests pass: 20/20 custom, 58/58 verification.

**Files**:
- `fds_c_interface.f90`: gate already removed

## Dependency Graph

```
[1. GET_LINKED_FV]     ← ✅ DONE (predictor/corrector barrier lambdas)
[2. CC_NO_FLUX_TS]     ← ✅ DONE (baroclinic_kernel_task.h)
       |
       v
[3. CC_NO_FLUX solve]  ← ✅ DONE (pressure_iteration_tasks.h)
[4. CC_VEL_ERROR_TS]   ← ✅ DONE (velocity_error_task.h)
[5. CC_EXCHANGE_TS]    ← ✅ DONE (fds_cc_exchange_prepare_fn in baroclinic_kernel_task.h)
       |
       v
[6. Enable subgraph]   ← ✅ DONE (gate removed, all items complete)
       |
       v
[7. Test]              ← ✅ 20/20 custom, 58/58 verification (tol=1e-6)
```

All items complete.

## Testing

CC_IBM test cases currently available:
- `shunn3_32_cc` — single mesh (won't trigger pressure subgraph)
- `sphere_helium_1mesh_cc` — single mesh
- `two_spheres_cc` — single mesh
- `sphere_helium_3meshes_cc` — 3 meshes, UGLMAT (can exercise pressure subgraph)

## Risk Assessment

| Item | Risk | Status |
|------|------|--------|
| CC_NO_FLUX thread safety | Low | ✅ Done — integrated, tests pass |
| CC_COMPUTE_VELOCITY_ERROR thread safety | Low | ✅ Done — integrated, tests pass |
| CC exchange FN_OMESH prep | Low | ✅ Done — kernel-side prep, no exchange modification |
| Multi-mesh CC_IBM correctness | Medium | ✅ sphere_helium_3meshes_cc passes (20/20, 58/58) |
| MPI CC exchange | High | Out of scope; same-rank only initially |

## Notes

- The existing CC_IBM multi-mesh segfault (mentioned in MEMORY.md) may be
  unrelated to the pressure subgraph — it could be in the predictor/corrector
  CC_IBM paths.  Must be fixed independently before this work can be tested.
- `CC_MATVEC_DEFINED` check in `MESH_CC_EXCHANGE` — must verify this flag is
  true when cut-cells are present with IBM forcing enabled.
- `GET_PRES_CFACE_BCS(NM, T, DT)` is called inside `ULMAT_SOLVER` (pres.f90:886),
  which is already in the solve kernel.  No separate integration needed.
- `GET_H_GUARD_CUTCELL` and `GET_H_CUTFACES` are initialization-time calls
  (inside `ULMAT_SOLVER` setup), not per-iteration.  No integration needed.
