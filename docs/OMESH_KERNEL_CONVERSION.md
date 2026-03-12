# OMESH Routine Kernel Conversion Plan

## Background

Many routines were kept sequential in orchestrators because they "read OMESH". Analysis shows OMESH is part of MESH_TYPE (`MESHES(NM)%OMESH(NOM)`) — it stores **per-mesh local copies** of neighboring mesh data, populated during MESH_EXCHANGE barriers. Two threads processing different meshes access completely separate memory. The actual obstacle was the `CALL POINT_TO_MESH(NM)` pattern (module-level mutable pointers), not OMESH itself.

### Data access patterns confirmed safe for parallel execution

| Access Pattern | Example | Race? |
|----------------|---------|-------|
| Read `M%OMESH(NOM)%<array>` | Ghost cell interpolation | No — NM's local copy |
| Write `M%OMESH(NOM)%<array>` | MATCH_VELOCITY averaging | No — NM's local copy |
| Read `MESHES(NOM)%DX/DY/DZ` | Neighbor grid geometry | No — read-only after init |
| Read `MESHES(NOM)%CELL_INDEX` | Neighbor topology | No — read-only after init |
| Write `M%<array>(II,JJ,KK)` | Own ghost cells | No — only NM's data |

## Conversion Method

For each candidate routine:

1. **Parallelizability audit** — Confirm POINT_TO_MESH is the *only* obstacle. Check for:
   - Writes to global arrays/accumulators (Q_DOT, M_DOT, Q_DOT_SUM, etc.)
   - Writes to other mesh data (`MESHES(NOM)%WALL(...)`, `MESHES(BACK_MESH)%...`)
   - Global counters (PARTICLE_TAG)
   - Cross-mesh particle/data buffer writes (PARTICLE_SEND_BUFFER)
   - RANDOM_NUMBER calls (compiler-dependent thread safety)
   - Subroutine calls that themselves have hidden global state

2. **Kernel conversion** — Replace POINT_TO_MESH usage:
   - Add `TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M` parameter (or POINTER for predictor/corrector phase selection)
   - Replace module pointers with `M%<component>` (e.g., `WALL` → `M%WALL`, `EXTERNAL_WALL` → `M%EXTERNAL_WALL`)
   - Replace `OMESH(NOM)` with `M%OMESH(NOM)`
   - Replace `MESHES(NOM)%DX` etc. with direct access (already safe — read-only geometry)
   - Remove `CALL POINT_TO_MESH(NM)`

3. **Integration** — Move from orchestrator sequential loop into parallel kernel task, or merge into existing kernel

4. **Verification** — Byte-identical test on 1–5 mesh configurations

---

## Group 1: Simple OMESH Ghost-Cell Fill Routines

These follow an identical pattern: loop over `N_EXTERNAL_WALL_CELLS`, read from `OMESH(NOM)%<array>`, average, write to local ghost cell. All are small, self-contained, and straightforward to convert.

### 1.1 VISCOSITY_BC

- **File**: `velo.f90:52–110` (58 lines)
- **Currently sequential in**: DivSetup orchestrators (pred + corr)
- **Reads**: `OMESH(NOM)%MU`, `OMESH(NOM)%KRES`, `OMESH(NOM)%D/DS`
- **Writes**: own `MU(II,JJ,KK)`, `KRES(II,JJ,KK)`, `D/DS(II,JJ,KK)` ghost cells
- **Other accesses**: `N_EXTERNAL_WALL_CELLS`, `WALL`, `EXTERNAL_WALL`, `BOUNDARY_COORD` (all module pointers from POINT_TO_MESH)

#### Parallelizability audit
- [x] No writes to global arrays
- [x] No writes to other mesh data
- [x] No global counters
- [x] No cross-mesh buffer writes
- [x] No RANDOM_NUMBER calls
- [x] No subroutine calls with hidden global state
- [x] **Verdict**: PARALLELIZABLE — pure per-mesh ghost-cell fill

#### Status: COMPLETE — Kernel in `velo_kernels.f90:VISCOSITY_BC_KERNEL`, integrated into DivSetup kernel task, verified byte-identical (1–5 mesh)

---

### 1.2 COMBUSTION_BC

- **File**: `fire.f90:808–843` (35 lines)
- **Currently sequential in**: CorrDivPart1 orchestrator
- **Reads**: `OMESH(NOM)%Q`
- **Writes**: own `Q(BC%II,BC%JJ,BC%KK)` ghost cells
- **Other accesses**: `N_EXTERNAL_WALL_CELLS`, `WALL`, `EXTERNAL_WALL`, `BOUNDARY_COORD`

#### Parallelizability audit
- [x] No writes to global arrays
- [x] No writes to other mesh data
- [x] No global counters
- [x] No cross-mesh buffer writes
- [x] No RANDOM_NUMBER calls
- [x] No subroutine calls with hidden global state
- [x] **Verdict**: PARALLELIZABLE — pure per-mesh ghost-cell fill

#### Status: COMPLETE — Kernel in `fire_kernels.f90:COMBUSTION_BC_KERNEL`, integrated into CorrDivPart1 kernel task, verified byte-identical (1–5 mesh)

---

### 1.3 ASSIGN_GHOST_VALUE

- **File**: `wall.f90:170–295` (125 lines)
- **Currently sequential in**: WallBC orchestrator (via WALL_BC_PREPROCESSING) and PredWallDiv orchestrator (via monolithic WALL_BC)
- **Reads**: `OMESH(NOM)%RHO/RHOS`, `OMESH(NOM)%ZZ/ZZS`, `MESHES(NOM)%DX/DY/DZ` (read-only geometry), `MESHES(NOM)%CELL_INDEX` (read-only topology)
- **Writes**: own `RHOP(II,JJ,KK)`, `ZZP(II,JJ,KK)`, `TMP(II,JJ,KK)`, `RSUM(II,JJ,KK)` ghost cells
- **Other accesses**: `EXTERNAL_WALL`, `WALL`, `BOUNDARY_COORD`, `BOUNDARY_PROP1`, `CELL_INDEX`, `CELL` — all module pointers
- **Note**: Called from within WALL_BC_PREPROCESSING loop, not standalone. Contains a section (line ~241–295) that reads `MESHES(NOM)%CELL_INDEX` and `MESHES(NOM)%CELL` for second-order interpolated boundary — these are read-only topology.

#### Parallelizability audit
- [x] No writes to global arrays
- [x] No writes to other mesh data (confirmed: MESHES(NOM)%DX/DY/DZ, CELL_INDEX, CELL are read-only geometry/topology)
- [x] No global counters
- [x] No cross-mesh buffer writes
- [x] No RANDOM_NUMBER calls
- [x] No subroutine calls with hidden global state (GET_SPECIFIC_GAS_CONSTANT is pure)
- [x] **Verdict**: PARALLELIZABLE — per-mesh ghost-cell interpolation with read-only neighbor geometry access

#### Status: COMPLETE — Kernel in `wall_kernels.f90:ASSIGN_GHOST_VALUE_KERNEL`, integrated into WallBC preprocessing kernel (C wrapper `C_FDS_WALL_BC_PREPROCESSING_KERNEL`), verified byte-identical (1–5 mesh)

---

## Group 2: Larger OMESH Routines

### 2.1 MATCH_VELOCITY

- **File**: `velo.f90:818–1057` (240 lines)
- **Currently sequential in**: PredFinal + CorrFinal orchestrators
- **Reads**: `OMESH(NOM)%U/US/V/VS/W/WS`, `MESHES(NOM)%DX/DY/DZ` (read-only geometry)
- **Writes**: own `UU/VV/WW` boundary cells, own `OMESH(NOM)%US/VS/WS` (NM's local copy — averaging at line 928), own `UVW_SAVE`, `U_GHOST/V_GHOST/W_GHOST`
- **CC_IBM path**: reads `MESHES(NOM)%FCVAR`, `MESHES(NOM)%CUT_FACE` (read-only geometry)
- **CC_IBM dispatch**: calls `CC_MATCH_VELOCITY(NM, PREDICTOR, .TRUE.)` — handled in C wrapper (dispatches to CC_MATCH_VELOCITY when CC_IBM, otherwise calls kernel)

#### Parallelizability audit
- [x] No writes to global arrays (T_USED timing skipped in kernel)
- [x] OMESH(NOM) writes are to NM's own local copy (not the actual neighbor)
- [x] MESHES(NOM) reads are all read-only geometry (DX, DY, DZ, FCVAR, CUT_FACE)
- [x] No global counters
- [x] No RANDOM_NUMBER calls
- [x] CC_IBM path handled at C wrapper level (CC_MATCH_VELOCITY still called sequentially when CC_IBM active; non-CC_IBM path uses kernel)
- [x] **Verdict**: PARALLELIZABLE — per-mesh velocity matching with read-only neighbor geometry and own-copy OMESH writes

#### Status: COMPLETE — Kernel in `velo_kernels.f90:MATCH_VELOCITY_KERNEL`, C wrapper `fds_match_velocity_kernel` handles CC_IBM dispatch, integrated into VelocityBCEdgesTask, verified byte-identical (1–5 mesh)

---

### 2.2 VELOCITY_BC_PREPROCESSING

- **File**: `velo.f90:710–784` (75 lines)
- **Currently sequential in**: PredFinal + CorrFinal orchestrators
- **Reads**: `M%OMESH(EWC%NOM)%US/VS/WS/U/V/W` (already uses `M%OMESH` directly!)
- **Writes**: own `UU/VV/WW` boundary cells, `M%DRAG_UVWMAX = 0`
- **Other accesses**: `N_EXTERNAL_WALL_CELLS`, `EXTERNAL_WALL` (were module pointers from POINT_TO_MESH)
- **Note**: Was *almost* already a kernel. Only needed POINT_TO_MESH removal + 2 module pointer replacements.

#### Parallelizability audit
- [x] No writes to global arrays
- [x] No writes to other mesh data
- [x] No global counters
- [x] No cross-mesh buffer writes
- [x] No RANDOM_NUMBER calls
- [x] No subroutine calls with hidden global state
- [x] **Verdict**: PARALLELIZABLE — pure per-mesh ghost-cell velocity fill

#### Status: COMPLETE — Modified in-place in `velo.f90:VELOCITY_BC_PREPROCESSING` (removed POINT_TO_MESH, replaced N_EXTERNAL_WALL_CELLS/EXTERNAL_WALL with M%), integrated into VelocityBCEdgesTask, verified byte-identical (1–5 mesh)

---

### 2.3 WALL_BC_PREPROCESSING

- **File**: `wall.f90:1502–1552` (50 lines)
- **Currently sequential in**: WallBC orchestrator
- **Calls**: ASSIGN_GHOST_VALUE (Group 1.3), NEAR_SURFACE_GAS_VARIABLES_KERNEL (already a kernel), HEAT_TRANSFER_COEFFICIENT (already thread-safe, converted to M-based)
- **Global state**: `fds_compute_wall_bc_dt_bc` and `fds_check_call_ht_1d` compute DT_BC and CALL_HT_1D from global `BC_CLOCK`/`WALL_COUNTER` — these single computations must remain sequential, but the per-mesh loop can be parallel.

#### Parallelizability audit
- [x] Global DT_BC/CALL_HT_1D computation isolated (already done in orchestrator)
- [x] ASSIGN_GHOST_VALUE is parallelizable (see 1.3 — ASSIGN_GHOST_VALUE_KERNEL complete)
- [x] NEAR_SURFACE_GAS_VARIABLES_KERNEL already thread-safe
- [x] HEAT_TRANSFER_COEFFICIENT already thread-safe (M-based)
- [x] No other global state writes
- [x] **Verdict**: PARALLELIZABLE — all sub-routines are per-mesh kernels

#### Status: COMPLETE — Full preprocessing loop in `C_FDS_WALL_BC_PREPROCESSING_KERNEL` (fds_c_interface.f90), integrated into WallBC kernel task, verified byte-identical (1–5 mesh)

---

## Group 3: Routines Without OMESH (POINT_TO_MESH Only)

### 3.1 SET_BAROCLINIC_FALSE

- **File**: `fds_c_interface.f90:642–645` (one-liner)
- **Currently sequential in**: DivSetup orchestrators (pred + corr)
- **Code**: `MESHES(NM)%BAROCLINIC_TERMS_ATTACHED = .FALSE.`

#### Parallelizability audit
- [x] Trivially per-mesh — writes only to MESHES(NM) component
- [x] **Verdict**: PARALLELIZABLE (trivial)

#### Status: COMPLETE — Already thread-safe (writes directly to MESHES(NM)), moved from DivSetup orchestrators to DivSetupKernelTask, verified byte-identical (1–5 mesh)

---

### 3.2 AGGLOMERATION / CALC_AGGLOMERATION

- **File**: `soot.f90:265–391` (126 lines)
- **Currently sequential in**: CorrDivSetup orchestrator (only if `N_AGGLOMERATION_SPECIES > 0`)
- **Pattern**: POINT_TO_MESH + pure per-mesh 3D loop over IBAR,JBAR,KBAR
- **Reads**: module pointers (U, V, W, RHO, TMP, ZZ, etc.)
- **Writes**: module pointers (ZZ — species mass fractions)
- **No OMESH access**

#### Parallelizability audit
- [x] No writes to global arrays (only writes to ZZ per-mesh species mass fractions)
- [x] No OMESH or cross-mesh access
- [x] No global counters
- [x] No RANDOM_NUMBER calls
- [x] GET_VISCOSITY: pure function (reads ZZ_GET, TMP), CUNNINGHAM: pure function
- [x] **Verdict**: PARALLELIZABLE — pure per-mesh 3D loop

#### Status: COMPLETE — Modified in-place in `soot.f90:CALC_AGGLOMERATION` (added M parameter, removed POINT_TO_MESH, local pointer aliases shadow module pointers), moved from CorrDivSetup orchestrator to DivSetupKernelTask (corrector path only via `work->estimated` check), verified byte-identical (1–5 mesh)

---

### 3.3 SYNTHETIC_TURBULENCE

- **File**: `turb.f90:723–930` (207 lines)
- **Currently sequential in**: PredFinal orchestrator (only if `SYNTHETIC_EDDY_METHOD`)
- **Reads**: `VENTS` (per-mesh, module pointer), `XC`, `YC` (module pointers)
- **Writes**: `VT%U_EDDY`, `VT%V_EDDY`, `VT%W_EDDY` (per-mesh vent data), `VT%X_EDDY`, `VT%Y_EDDY`, `VT%Z_EDDY` (eddy positions)
- **Calls**: `EDDY_POSITION` → `RANDOM_NUMBER`, `EDDY_AMPLITUDE` → `RANDOM_NUMBER`

#### Parallelizability audit
- [x] No writes to global arrays
- [x] Writes are all to per-mesh VENTS data
- [x] No cross-mesh access
- [ ] **RANDOM_NUMBER thread safety** — Fortran intrinsic RANDOM_NUMBER uses compiler-dependent RNG state. gfortran uses per-thread seed (thread-safe), ifort uses global seed (NOT thread-safe). Parallelizing would change random number sequence, breaking byte-identical reproducibility.
- [x] EVALUATE_RAMP: reads RAMPS (global read-only array) — safe
- [x] **Verdict**: CONDITIONALLY PARALLELIZABLE — safe for data access, but RANDOM_NUMBER thread safety is compiler-dependent. Cannot verify byte-identical since test cases don't use SYNTHETIC_EDDY_METHOD.

#### Status: DEFERRED — Kept sequential in PredFinal orchestrator. Would need compiler-specific verification and potentially per-mesh RNG seeding for deterministic behavior.

---

## Group 4: Routines Requiring Deep Analysis

### 4.1 CC_VELOCITY_BC

- **File**: `ccib_velocity.f90:2609–3699` (~1090 lines)
- **Currently sequential in**: DivSetup orchestrators (CC_IBM only), PredFinal/CorrFinal collectors (CC_IBM only)
- **Pattern**: Large CC_IBM routine, uses POINT_TO_MESH extensively, accesses `CUT_FACE`, `CC_EDGE`, many module-level pointers
- **Module pointers used**: US, VS, WS, U, V, W, ZZS, ZZ, RHOS, RHO, CUT_FACE, CUT_EDGE, CC_RCEDGE, CC_IBEDGE, EDGE, CELL, MU, TMP, DX, DY, DZ, IBAR, JBAR, KBAR, IBP1, JBP1, KBP1, FCVAR, DRAG_UVWMAX
- **CONTAINS subroutines** (4): CC_CUTEDGE_DUIDXJ_TAU_OMG, CC_RCEDGE_DUIDXJ, CC_RCEDGE_TAU_OMG, CC_EDGE_TAU_OMG

#### Parallelizability audit
- [x] All data accesses (CUT_FACE, CC_EDGE, CC_RCEDGE, CC_IBEDGE, FCVAR) are per-mesh via POINT_TO_MESH
- [x] No OMESH access (no cross-mesh data reads)
- [x] No writes to MESHES(NOM) or other mesh data
- [x] Subroutine calls (WALL_MODEL, GET_VISCOSITY) are pure/thread-safe
- [x] No global counters or accumulators
- [x] **Verdict**: PARALLELIZABLE after kernel conversion (POINT_TO_MESH is the only obstacle)

#### Status: DEFERRED — Parallelizable in principle, but requires ~1090 line kernel conversion using Pattern 3 (local pointer alias shadowing). For non-CC_IBM runs, the call is a complete no-op (guarded by `IF (CC_IBM)`). DivSetup orchestrators and PredFinal/CorrFinal collectors conditionally eliminated for non-CC_IBM via `fds_is_cc_ibm()` check at graph construction time.

---

### 4.2 UPDATE_GLOBAL_OUTPUTS

- **File**: `dump.f90:64–83` (20 lines wrapper)
- **Currently sequential in**: CorrFinal collector
- **Calls**:
  - `UPDATE_HRR(DT,NM)`: writes to `Q_DOT` (global array), `M_DOT` (global array), `Q_DOT_SUM`/`M_DOT_SUM` (global accumulators), `ENTHALPY_SUM(NM)` (per-mesh slot)
  - `UPDATE_MASS(DT,NM)`: similar global accumulator pattern
  - `UPDATE_FIRE_SPREAD_OUTPUTS(T,DT,NM)`: per-mesh fire spread data
  - `UPDATE_DEVICES_1(T,DT,NM)`: device state updates (~530 lines), writes to `DEVICE%PRIOR_STATE` (global DEVC array)

#### Parallelizability audit
- [x] `Q_DOT`, `M_DOT` — global arrays written by all meshes → **race condition**
- [x] `Q_DOT_SUM`, `M_DOT_SUM` — global accumulators → **race condition**
- [x] `ENTHALPY_SUM(NM)` — per-mesh, safe
- [x] `UPDATE_DEVICES_1` — writes to `DEVICE%PRIOR_STATE`, `DEVICE%SMOOTHED_VALUE` (global DEVC array) → **race condition**
- [x] **Verdict**: MUST REMAIN SEQUENTIAL — global accumulators and device state mutations

#### Status: MUST REMAIN SEQUENTIAL — Kept in CorrFinal collector. Cannot be parallelized due to global Q_DOT/M_DOT accumulators and DEVICE state writes. A partial-sums approach (like CorrRadiation RAD_Q_SUM) is theoretically possible for Q_DOT/M_DOT but not for DEVICE state.

---

### 4.3 CC_MATCH_VELOCITY

- **File**: `ccib_velocity.f90:1921–2500` (~580 lines)
- **Currently called by**: `fds_match_velocity_kernel` C wrapper (dispatches to CC_MATCH_VELOCITY when CC_IBM is active, otherwise calls MATCH_VELOCITY_KERNEL)
- **Pattern**: CC_IBM version of MATCH_VELOCITY, uses POINT_TO_MESH, reads OMESH, accesses CUT_FACE

#### Parallelizability audit
- [x] OMESH access patterns: reads `M%OMESH(NOM)%` data (per-mesh local copy, safe)
- [x] CUT_FACE writes: writes to `MESHES(NOM)%CUT_FACE(ICF)%VELS_OMESH` and `%VEL_OMESH` and `%VEL_LNK_OMESH` — **cross-mesh writes to neighbor's CUT_FACE data**
- [x] **Verdict**: MUST REMAIN SEQUENTIAL — cross-mesh CUT_FACE writes (two threads processing different meshes could write to the same neighbor's CUT_FACE simultaneously)

#### Status: MUST REMAIN SEQUENTIAL — Kept sequential in `fds_match_velocity_kernel` C wrapper (dispatched via POINT_TO_MESH for CC_IBM). Cross-mesh `MESHES(NOM)%CUT_FACE` writes prevent parallelization. The non-CC_IBM path (MATCH_VELOCITY_KERNEL) is already parallel.

---

## Group 5: Genuinely Sequential (No Conversion Planned)

These routines have confirmed cross-mesh writes or global state mutations that prevent parallelization.

### 5.1 WALL_BC_FINALIZE

- **File**: `wall.f90:1432–1489` (57 lines)
- **Reason**: HAS_BACK_MESH cells read/write `MESHES(BACK_MESH)%WALL/BOUNDARY_PROP1/BOUNDARY_PROP2` — genuine cross-mesh writes. DEPOSIT_PARTICLE_MASS writes to global `Q_DOT`, `M_DOT`.
- **Status**: MUST REMAIN SEQUENTIAL

### 5.2 INSERT_ALL_PARTICLES

- **File**: `part.f90:119–1801` (1682 lines)
- **Reason**: Global `PARTICLE_TAG` counter (incremented by NMESHES per particle). `OMESH%PARTICLE_SEND_BUFFER` writes for cross-mesh particle injection.
- **Status**: MUST REMAIN SEQUENTIAL

### 5.3 MOVE_PARTICLES

- **File**: `part.f90:1806–3455` (1649 lines)
- **Reason**: `ADD_TO_PARTICLE_SEND_BUFFER` for particles leaving mesh boundaries — writes to `OMESH%PARTICLE_SEND_BUFFER` which is consumed by MESH_EXCHANGE.
- **Status**: MUST REMAIN SEQUENTIAL

### 5.4 PARTICLE_MASS_ENERGY_TRANSFER

- **File**: `part.f90:3463–4553` (1090 lines)
- **Currently sequential in**: CorrParticle orchestrator
- **Needs audit**: May be per-mesh only (no particle buffer writes). If so, could be promoted to Group 3.
- **Status**: NEEDS AUDIT (deferred — low priority, particle routines are fast)

---

## Special Case: Predictor WALL_BC Decomposition

### Current State

In the **corrector** phase, WALL_BC is already decomposed into three phases:
- `WALL_BC_PREPROCESSING` (sequential in orchestrator)
- `WALL_BC_PROCESS_CELLS_KERNEL` (parallel kernel task)
- `WALL_BC_FINALIZE` (sequential in collector)

In the **predictor** phase (`PredWallDiv`), the monolithic `WALL_BC(T,DT,NM)` is called, which internally calls all three phases sequentially. This should be decomposed to match the corrector pattern.

### Action Item

- [x] Replace `fds_wall_bc(md->t, md->dt, md->nm)` in PredWallDivOrchestrator with the three-phase decomposition
- [x] Share the same kernel task and sub-graph pattern as the corrector WallBC

#### Status: COMPLETE — Reused corrector's WallBC sub-graph (`buildWallBCSubgraph`) in predictor pipeline. Removed PredWallDivOrchestrator (was only calling monolithic WALL_BC). `fds_check_call_ht_1d` correctly returns 0 during predictor (checks CORRECTOR flag). Verified byte-identical (1–5 mesh)

---

## Integration Strategy

### Phase 1: Low-hanging fruit (Groups 1 + 3.1) — ✅ COMPLETE (Group 1)

Group 1 complete: VISCOSITY_BC_KERNEL, COMBUSTION_BC_KERNEL, ASSIGN_GHOST_VALUE_KERNEL all converted and verified byte-identical on 1–5 mesh tests. Sequential preprocessing removed from DivSetup (pred+corr), CorrDivPart1, and WallBC orchestrators. SET_BAROCLINIC_FALSE (Group 3.1) remains.

**Impact achieved**: DivSetup orchestrators no longer run VISCOSITY_BC sequentially. CorrDivPart1 orchestrator no longer runs COMBUSTION_BC sequentially. WallBC orchestrator no longer runs WALL_BC_PREPROCESSING sequentially — all moved to parallel kernel tasks.

### Phase 2: OMESH routines (Groups 1.3 + 2) — ✅ COMPLETE

All converted: ASSIGN_GHOST_VALUE_KERNEL, MATCH_VELOCITY_KERNEL, VELOCITY_BC_PREPROCESSING (in-place fix), WALL_BC_PREPROCESSING_KERNEL. All verified byte-identical on 1–5 mesh tests.

**Impact achieved**: CorrFinal orchestrator has NO sequential preprocessing — dispatches work immediately. PredFinal orchestrator only keeps SYNTHETIC_TURBULENCE (RANDOM_NUMBER concern, Group 3.3). WallBC orchestrator only computes global DT_BC/CALL_HT_1D then dispatches immediately.

### Phase 3: Groups 3 + deeper analysis — ✅ COMPLETE (Groups 3.1 + 3.2), DEFERRED (3.3)

SET_BAROCLINIC_FALSE and CALC_AGGLOMERATION moved to parallel DivSetupKernelTask. SYNTHETIC_TURBULENCE deferred (RANDOM_NUMBER thread safety concern, kept sequential). Remaining Group 4 items (CC_VELOCITY_BC, UPDATE_GLOBAL_OUTPUTS) need deeper analysis.

**Impact achieved**: DivSetup orchestrators (pred+corr) reduced to only CC_VELOCITY_BC sequential call (CC_IBM-specific, Group 4).

### Phase 4: Predictor WALL_BC decomposition — ✅ COMPLETE

Reused corrector's `buildWallBCSubgraph` in predictor pipeline. Removed PredWallDivOrchestrator (only existed to call monolithic `fds_wall_bc`). Predictor WALL_BC now runs ~90% of wall cells in parallel via WallBCKernelTask.

**Impact achieved**: Predictor pipeline no longer has any monolithic sequential WALL_BC call. Both predictor and corrector use the same three-phase WallBC sub-graph.

### Phase 5: Group 4 deep analysis + conditional graph construction — ✅ COMPLETE

**Audit results**:
- CC_VELOCITY_BC: PARALLELIZABLE in principle (POINT_TO_MESH only obstacle), but ~1090 line kernel conversion deferred. For non-CC_IBM, call is a complete no-op.
- UPDATE_GLOBAL_OUTPUTS: MUST REMAIN SEQUENTIAL (global Q_DOT/M_DOT accumulators, DEVICE state)
- CC_MATCH_VELOCITY: MUST REMAIN SEQUENTIAL (cross-mesh MESHES(NOM)%CUT_FACE writes)

**Conditional graph construction**: Added `fds_is_cc_ibm()` C interface to query CC_IBM flag at graph construction time. For non-CC_IBM runs, the following synchronization barriers are eliminated:
- DivSetup orchestrators (pred+corr): MeshExchange dispatches directly to parallel kernel task
- VelocityPredictor orchestrator+collector: PressureIteration dispatches directly to parallel kernel task
- VelocityCorrector orchestrator+collector: PressureIteration dispatches directly to parallel kernel task
- PredFinal collector: kernel task outputs directly (no CC_VELOCITY_BC post-processing)

**Data type cleanup**: Eliminated VelocityPredictorWork, VelocityCorrectorWork, VelocityBCWork intermediate types. All kernel tasks now use MeshData → MeshData directly. VelocityBCEdgesTask takes `applyToEstimated` as constructor parameter.

**Impact achieved**: For non-CC_IBM runs (standard FDS cases), 7 unnecessary synchronization barriers removed. Only inherently global barriers remain (MESH_EXCHANGE, PRESSURE_ITERATION, HVAC, COMBUSTION, INSERT_ALL_PARTICLES, MOVE_PARTICLES, WALL_BC_FINALIZE, UPDATE_GLOBAL_OUTPUTS, SYNTHETIC_TURBULENCE).

---

## Summary

| Group | Routines | Lines | Outcome |
|-------|----------|-------|---------|
| 1 (simple OMESH) | VISCOSITY_BC, COMBUSTION_BC, ASSIGN_GHOST_VALUE | ~218 | ✅ Merged into parallel kernels |
| 2 (larger OMESH) | MATCH_VELOCITY, VELOCITY_BC_PREPROCESSING, WALL_BC_PREPROCESSING | ~365 | ✅ Removed orchestrator barriers |
| 3 (POINT_TO_MESH only) | SET_BAROCLINIC_FALSE, AGGLOMERATION, SYNTHETIC_TURBULENCE | ~334 | ✅ Groups 3.1+3.2 parallel; 3.3 deferred (RANDOM_NUMBER) |
| 4 (deep analysis) | CC_VELOCITY_BC, UPDATE_GLOBAL_OUTPUTS, CC_MATCH_VELOCITY | ~1690 | ✅ Audited: CC_VELOCITY_BC deferred (parallelizable but large), others sequential. Non-CC_IBM barriers eliminated. |
| 5 (genuinely sequential) | WALL_BC_FINALIZE, INSERT_ALL_PARTICLES, MOVE_PARTICLES | ~3388 | No conversion (cross-mesh writes) |
