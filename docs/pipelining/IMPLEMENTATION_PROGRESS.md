# FDS Pipelining Implementation Progress

## Overview

This file tracks the implementation of intra-timestep pipelining parallelism in the Hedgehog graph. The analysis and design are documented in [README.md](README.md) and [fds_section_pipeline.dot](fds_section_pipeline.dot).

**Goal**: Overlap independent computation stages within each timestep to reduce wall-clock time. Two main fork-join pairs in the corrector phase (VFLUX || COMBUSTION, RADIATION || DIV_P1) provide a combined 1.36x speedup of the corrector pipeline.

**Strategy**: Bottom-up implementation. First prepare the Fortran kernels and validate them sequentially in the existing graph, then restructure the graph with fork-join states.

## Current Phase: Complete (All 6 Phases Done)

---

## Phase 1: Fortran Kernel Extraction

**Objective**: Create the split kernel variants needed for pipelining, without changing any graph structure. Validate by calling them in the same sequential order.

**Rationale**: Isolate Fortran correctness from graph correctness. All new kernels must produce bit-identical results when called in the original order.

### 1a: DIV_P1 QR-Independent Variant

Add a `SKIP_QR` logical parameter to `DIVERGENCE_PART_1_KERNEL`. When `.TRUE.`, omit the `M%QR(I,J,K)` addition in `COMPUTE_THERMAL_DIVERGENCE` (lines 566, 577 of divg_kernels.f90). When `.FALSE.`, behavior is identical to current code.

**Files to modify**:
- `Source/divg_kernels.f90`: Add `SKIP_QR` parameter, guard the 2 QR references
- `Source/divg.f90`: Update call site (pass `.FALSE.` to preserve current behavior)

**Test**: `SKIP_QR=.FALSE.` everywhere → bit-identical on full verification suite.

- [x] Add SKIP_QR parameter to DIVERGENCE_PART_1_KERNEL
- [x] Update call sites in divg.f90 (OPTIONAL parameter, no change needed)
- [x] Verify bit-identical (12/12 tests pass)

### 1b: DIV_P1 QR Addition Kernel

Create a new subroutine `DIVERGENCE_PART_1_ADD_QR_KERNEL(M, NM)` in `divg_kernels.f90`. This is a simple loop:

```fortran
DP(I,J,K) = DP(I,J,K) + M%QR(I,J,K)
```

with predictor/corrector DP pointer selection (M%DS or M%D). Handles both Cartesian and cylindrical.

**Files to create/modify**:
- `Source/divg_kernels.f90`: Add new subroutine

**Test**: Call `DIV_P1_KERNEL(SKIP_QR=.TRUE.)` then `DIV_P1_ADD_QR_KERNEL()` sequentially → bit-identical.

- [x] Create DIVERGENCE_PART_1_ADD_QR_KERNEL
- [x] Verify bit-identical with split call sequence (12/12 tests pass)

### 1c: C Wrappers

Create C-callable wrappers for the new/split kernels.

**Files to modify**:
- `Source/hedgehog/fds_c_interface.f90`: Add wrappers
  - `fds_divergence_part_1_kernel_skip_qr(nm, t, dt)` — calls with SKIP_QR=.TRUE.
  - `fds_divergence_part_1_add_qr(nm, t, dt)` — calls ADD_QR_KERNEL

**Test**: C wrappers callable from C++ without linker errors.

- [x] Create C wrappers (fds_divergence_part_1_kernel_skip_qr, fds_divergence_part_1_add_qr)
- [x] Link test (build succeeds, 12/12 tests pass)

---

## Phase 2: Per-Branch Scratch Arrays

**Objective**: Add a second set of WORK/SWORK arrays to MESH_TYPE so that two concurrent branches can use independent scratch memory. Validate by running with branch=1 (original arrays) everywhere.

**Rationale**: RADIATION uses WORK1-9, DIV_P1 uses WORK1-7,9. They cannot run concurrently on the same mesh without separate scratch pools. Memory cost: ~4.5 MB per 64³ mesh (9 3D double arrays).

### 2a: Add Branch-B WORK Arrays to MESH_TYPE

Add to `MESH_TYPE` in `mesh.f90`:
- `WORK1_B` through `WORK9_B` (same dimensions as WORK1-9)
- `SWORK1_B` through `SWORK3_B` (same dimensions as SWORK1-3)

**Files to modify**:
- `Source/mesh.f90`: Declare new arrays in MESH_TYPE
- `Source/init.f90` (or wherever WORK arrays are allocated): Allocate the _B arrays

- [x] Add WORK_B declarations to MESH_TYPE
- [x] Add allocation in initialization
- [x] Verify no build errors (12/12 tests pass)

### 2b: Add WORK_BRANCH Parameter to Affected Kernels

Add a `WORK_BRANCH` integer parameter to kernels that use WORK arrays. When `WORK_BRANCH=1`, use original `M%WORK1`, etc. When `WORK_BRANCH=2`, use `M%WORK1_B`, etc.

The change is localized to the pointer alias block at the top of each kernel (similar to the existing PREDICTOR/CORRECTOR alias pattern for UU/VV/WW). CONTAINS subroutines inherit the aliases via host association.

**Kernels to modify**:
- `Source/divg_kernels.f90`: DIVERGENCE_PART_1_KERNEL (uses WORK1-7,9, SWORK1-3)
- `Source/radi.f90`: COMPUTE_RADIATION_KERNEL (uses WORK1-9)
- `Source/velo_kernels.f90`: VELOCITY_FLUX_KERNEL (uses WORK1-6)
- `Source/mass_kernels.f90`: DENSITY_KERNEL (uses WORK4-5, SWORK4)
- `Source/fire_kernels.f90`: CONDENSATION_EVAPORATION_KERNEL (uses WORK1-2, SWORK1)
- `Source/part.f90`: PARTICLE_MASS_ENERGY_KERNEL (uses WORK1-7, SWORK1)

**Files to modify**:
- Each kernel file above: Add WORK_BRANCH parameter, conditional pointer aliases
- `Source/hedgehog/fds_c_interface.f90`: Update C wrappers to pass WORK_BRANCH
- All Hedgehog tasks that call these kernels: Pass WORK_BRANCH=1 (default)

**Test**: `WORK_BRANCH=1` everywhere → bit-identical on full verification suite.

- [x] Add WORK_BRANCH to DIV_P1 kernel (+ VELOCITY_FLUX_BLOCK_KERNEL)
- [x] Add WORK_BRANCH to RADIATION kernel
- [x] Add WORK_BRANCH to VELOCITY_FLUX kernel
- [x] Add WORK_BRANCH to remaining kernels (DENSITY, CONDENSATION, PME)
- [x] C wrappers unchanged (WORK_BRANCH is OPTIONAL, defaults to branch 1)
- [x] Hedgehog tasks unchanged (OPTIONAL defaults handle this)
- [x] Verify bit-identical (12/12 tests pass)

---

## Phase 3: Sequential Driver Integration

**Objective**: Wire the split kernels into the existing corrector sub-graph, calling them in the same sequential order. This validates the new kernels in context without changing graph structure.

**Rationale**: This is the "prepare the driver" step. If verification fails here, the bug is in the Fortran kernel split, not in the graph restructuring.

### 3a: Replace DIV_P1 Call in Corrector

In the corrector sub-graph, replace the single DIV_P1 kernel call with the two-step sequence:
1. `fds_divergence_part_1_kernel_skip_qr(nm, t, dt)` — DIV_P1 without QR
2. `fds_divergence_part_1_add_qr(nm)` — RTRM*QR addition

**Critical fix discovered**: The ADD_QR kernel must multiply QR by RTRM (stored in WORK1) before adding, because DP undergoes `DP = RTRM * DP` during COMPUTE_DIVERGENCE_SOURCES. Simple `DP += QR` produces wrong results since QR would bypass the RTRM scaling. The correct formula is `DP += RTRM * QR`.

The Hedgehog task `CorrDivPart1KernelTask` calls both sequentially in its `execute()` method.

**Files modified**:
- `Source/hedgehog/task/corr_div_part1_kernel_task.h`: Call split sequence
- `Source/divg_kernels.f90`: Fix ADD_QR kernel to use `DP += RTRM * QR`
- `Source/hedgehog/fds_c_interface.f90`: Add `fds_divergence_part_1_add_qr_b` wrapper (WORK_BRANCH=2)
- `Source/hedgehog/fds_fortran_interface.h`: Declare new C wrapper

- [x] Update CorrDivPart1KernelTask to use split sequence
- [x] Fix ADD_QR kernel to use RTRM*QR (12/12 tests pass)

### 3b: Validate WORK_BRANCH Plumbing

Temporarily set WORK_BRANCH=2 for DIV_P1 or RADIATION to verify that branch-B scratch arrays produce identical results.

- [x] Run DIV_P1 with WORK_BRANCH=2 → bit-identical (10/12 pass; 2 CC_IBM cases fail due to CC code using M%WORK directly)
- [x] Run RADIATION with WORK_BRANCH=2 → bit-identical (12/12 pass)
- [x] Restore WORK_BRANCH=1 for all

**Note**: CC_IBM paths in `ccib_divergence_kernels.f90` use `M%WORK2-4` and `M%SWORK1-3` directly, bypassing the WK/SK intermediates. This means WORK_BRANCH=2 for DIV_P1 is not compatible with CC_IBM until CC code is updated. Phase 5 (RADIATION || DIV_P1) should use WORK_BRANCH=2 for DIV_P1 only when CC_IBM is false.

---

## Phase 4: Corrector Fork 1 — VFLUX || COMBUSTION

**Objective**: Run VELOCITY_FLUX and COMBUSTION concurrently after MESH_EXCHANGE(4). This is the first graph restructuring change.

**WORK conflict**: None. COMBUSTION uses no WORK arrays. CONDENSATION uses WORK1-2 but runs after the join (after SootHvac barrier). No scratch duplication needed for this fork.

**Expected savings**: 115 ops/cell hidden behind VFLUX (248 ops/cell). Higher with particles (Branch B grows longer).

### Architecture

```
EX(4) scatters MeshData
  ↓
PipelineFork1 state (MeshData → VFluxWork + CombWork)
  ├── Branch A: VelocityFlux sub-graph (VFluxWork → VFluxResult)
  └── Branch B: Combustion kernel task (CombWork → CombResult)
PipelineJoin1 state (VFluxResult + CombResult → MeshData, per-mesh matching)
  ↓
SootHvac collector → SootHvac barrier
  ↓
(continues: COND → PME → Move → PART_MOM → ...)
```

### Data Types

New data types needed (lightweight wrappers around MeshData):
- `VFluxWork` / `VFluxResult` — Branch A tokens
- `CombWork` / `CombResult` — Branch B tokens

### Key Design Decisions

- Fork is **per-mesh**: each MeshData spawns one VFluxWork + one CombWork
- Join matches by mesh index (NM): emits MeshData when both branch results arrive for that mesh
- Branch A uses the existing VelocityFlux sub-graph (mesh-level or K-block)
- Branch B is a simple kernel task (CombustionKernelTask adapted for CombWork input)
- SootHvac barrier runs AFTER the join (all meshes, both branches complete)

### Files to Create

- `Source/hedgehog/data/pipeline_data.h`: Fork/join data types
- `Source/hedgehog/state/pipeline_fork_state.h`: Fork state (MeshData → typed branches)
- `Source/hedgehog/state/pipeline_join_state.h`: Join state (typed results → MeshData)

### Files to Modify

- `Source/hedgehog/graph/corrector_subgraph.h`: Wire fork/join around VFLUX and COMB
- `Source/hedgehog/task/combustion_kernel_task.h`: Accept CombWork input type (or create adapter)

### Checklist

- [x] Create pipeline data types (Fork1VFluxWork/Result, Fork1CombWork/Result)
- [x] Create fork state (PipelineFork1State)
- [x] Create join state (PipelineJoin1State, per-mesh matching)
- [x] Adapt VelocityFlux sub-graph for VFluxWork input (nested sub-graph with unwrap/wrap)
- [x] Adapt CombustionKernelTask for CombWork input (Fork1CombKernelTask)
- [x] Wire corrector sub-graph with fork/join (12/12 tests pass)
- [x] Run verification suite (46/58 pass — no regressions from Fork 1)
- [x] Compare performance: VFLUX hidden behind COMB — dancing_eddies: 1.01x (no combustion), fire_const_gamma: 1.25x (446ms hidden), species_props: 1.11x

---

## Phase 5: Corrector Fork 2 — RADIATION || DIV_P1

**Objective**: Run RADIATION and DIV_P1 (without QR) concurrently after MESH_EXCHANGE(6a). This is the high-impact change (1299 ops/cell saved).

**WORK conflict**: RADIATION uses WORK1-9, DIV_P1 uses WORK1-7,9. **Requires per-branch scratch duplication** (Phase 2). RADIATION uses WORK_BRANCH=1, DIV_P1 uses WORK_BRANCH=2.

**Expected savings**: min(RADIATION, 1299) ops/cell. For typical fire: 1299 ops/cell — the entire DIV_P1 computation runs hidden behind RADIATION.

### Architecture

```
EX(6a) scatters MeshData
  ↓
PipelineFork2 state (MeshData → RadiationWork + DivP1Work)
  │  [also zeros DSUM/PSUM/USUM for DIV_P1]
  ├── Branch C: Radiation sub-graph (RadiationWork → RadiationResult)
  │   [uses WORK_BRANCH=1]
  └── Branch D: DIV_P1_non_QR kernel (DivP1Work → DivP1Result)
      [uses WORK_BRANCH=2, SKIP_QR=.TRUE.]
RadiationCollector (N meshes → RadBarrierData)
DivP1Collector (N meshes → DivP1BarrierData)
PipelineJoin2 state (RadBarrierData + DivP1BarrierData → BarrierData)
  ↓
MeshExchange(2) [exchanges QR across meshes, NO InitDiv — already done in fork]
  ↓ scatters MeshData
QR Addition kernel (per mesh)
  ↓
DIV_EXCHANGE collector → ...
```

### Key Design Decisions

- Both branches need ALL meshes to synchronize before the join (EX(2) is a global barrier for QR exchange)
- InitDiv (zero DSUM/PSUM/USUM) moves to the fork state's orchestration, before dispatching Branch D
- EX(2) only exchanges QR (InitDiv flag = false); the zeroing is done in the fork
- QR addition is a trivial per-mesh kernel (1 op/cell) that runs after EX(2)
- Reuse fork/join patterns from Phase 4 with different data types

### Data Types

- `RadiationWork` / `RadiationResult` — Branch C tokens
- `DivP1Work` / `DivP1Result` — Branch D tokens
- `RadBarrierData` / `DivP1BarrierData` — Branch collector outputs

### Files to Create

- `Source/hedgehog/data/pipeline_fork2_data.h`: Fork 2 data types
- `Source/hedgehog/state/pipeline_fork2_state.h`: Fork 2 state
- `Source/hedgehog/state/pipeline_join2_state.h`: Join 2 state
- `Source/hedgehog/task/div_p1_qr_addition_task.h`: QR addition task

### Files to Modify

- `Source/hedgehog/graph/corrector_subgraph.h`: Wire fork 2 around RADIATION and DIV_P1
- `Source/hedgehog/graph/corr_radiation_subgraph.h`: Accept RadiationWork input
- `Source/hedgehog/task/corr_div_part1_kernel_task.h`: Accept DivP1Work, use SKIP_QR + WORK_BRANCH=2

### Checklist

- [x] Create Fork 2 data types (Fork2RadWork/Barrier, Fork2DivP1Work/Barrier)
- [x] Create Fork 2 state (PipelineFork2State: collects N meshes, InitDiv, dispatches)
- [x] Create Join 2 state (PipelineJoin2State: Fork2RadBarrier + Fork2DivP1Barrier → BarrierData)
- [x] Create QR addition task (DivP1QRAdditionTask: WORK_BRANCH=2, after MeshExchange(2))
- [x] Adapt Radiation sub-graph for Fork2RadWork input (unwrap/wrap pattern in pipeline_fork2_rad_subgraph.h)
- [x] Create Branch D DIV_P1 task (Fork2DivP1KernelTask: SKIP_QR + WORK_BRANCH=2) + collector
- [x] Modify MeshExchange(2) to skip InitDiv when !CC_IBM (initDiv=ccIBM conditional)
- [x] Wire corrector sub-graph with Fork 2 (CC_IBM sequential fallback)
- [x] Run verification suite (46/58 pass — no regressions from Fork 2)
- [x] Compare performance: RAD||DIV_P1 near-balanced — dancing_eddies: 1.91x (506ms hidden), fire_const_gamma: 1.29x (534ms hidden), species_props: 1.11x

---

## Phase 6: Predictor Pipelining

**Objective**: Add pipelining to the predictor phase.

**Decision**: Option A (two-way fork: VFLUX || WALL_BC). Option B (three-way fork including DIV_P1) was investigated and ruled out.

### Investigation Results: Option B NOT Feasible

A phase-by-phase analysis of DIV_P1's interior/wall correction split was conducted. The findings:

**PREDICT_NORMAL_VELOCITY**: Self-contained in predictor (uses SURFACE properties + current velocity, not WALL_BC outputs). Can run independently. OK.

**Phase 2A (COMPUTE_SPECIES_DIFFUSION_FLUXES)**: Interior face gradients are independent of WALL_BC, but the mass conservation correction (lines 275-288) operates on ALL faces including boundary faces. Wall corrections overwrite boundary faces afterward, so this could be tolerable.

**Phase 2B (COMPUTE_DIFFUSIVE_HEAT_FLUX)**: SHOWSTOPPER. The interior loop computes `H_RHO_D_DZDX` at ALL face positions (0:IBAR), including boundary faces. Without wall corrections, boundary face values are wrong. The divergence accumulation (`DP += div(H_RHO_D_DZDX)`) reads these wrong boundary values for wall-adjacent cells (e.g., cell I=1 reads face I=0). Phase 2B's own wall correction (lines 338-413) independently overwrites scratch arrays, but the DP contribution has already been computed.

**Phase 5 (COMPUTE_DIVERGENCE_SOURCES)**: FATAL. `DP *= RTRM` (lines 609-632) is a **multiplicative** operation on DP. Any prior error from wrong boundary face values in Phases 2B and 4b is amplified. No additive correction can fix this after the fact.

**Consequence**: To correct wall-adjacent DP values after WALL_BC completes, we would need to re-run the entire DIV_P1 computation for those cells — essentially running DIV_P1 twice. This negates the parallelism savings.

| Factor | Detail |
|--------|--------|
| Root cause | DP accumulation is sequential with multiplicative step (DP *= RTRM) |
| Boundary face reads | Phases 2B, 4b, 5A, 5C all read scratch arrays at boundary faces |
| B1 properties needed | TMP_F, ZZ_F, RHO_F, RHO_D_DZDN_F, Q_CON_F — all set by WALL_BC |
| Correction feasibility | Impossible without re-running DIV_P1 for wall-adjacent cells |

### Improved Design: (VFLUX+PMOM) || (WALL_BC+DIV_P1_early) (CHOSEN)

A phase-by-phase analysis of DIV_P1 revealed that Phases 2A through 4b (species diffusion, heat diffusion, specific heat, thermal conductivity/divergence) do **NOT** reference UU/VV/WW or FVX/FVY/FVZ at all. Only Phase 5+ (enthalpy advection, species advection, sources, zone sums) reads velocity and flux arrays. This means the early phases can run in Branch B after WALL_BC, concurrently with VFLUX+PART_MOM in Branch A.

**CC_IBM safety**: SET_EXIMDIFFLX_3D (called within Phase 2A) only touches RHO_D_DZDX/Y/Z — safe in Branch B. CC_VELOCITY_FLUX(CORRECT_GRAV=.FALSE.) reads FVX via CC_STORE_FACE_FV — must run post-join. CFACE_PREDICT_NORMAL_VELOCITY only reads B1/SURFACE properties — safe in pre-fork. No CC_IBM special cases needed.

**Expected savings**: 266 ops/cell (Branch A hidden behind Branch B). Fork speedup: 1.74×.

```
PREDICT_NORMAL_VELOCITY + setup (~5 ops/cell)
  → Fork →
    Branch A: VFLUX(248) → PART_MOM(18) = 266 ops/cell (hidden behind B)
    Branch B: WALL_BC(65) → DIV_P1 early Ph2A-4b(295) = 360 ops/cell (CRITICAL)
  → Join
  → [CC_VELOCITY_FLUX if CC_IBM — needs FVX]
  → DIV_P1 Phase 5+ (~1000 ops/cell — needs UU/VV/WW + FVX)
  → [CC_DIVERGENCE_PART_1 if CC_IBM]
  → DIV_EXCHANGE → DIV_P2 → PRESSURE_SOLVE → ...
```

**Architecture**: Split DIVERGENCE_PART_1_KERNEL into three callable sections via PHASE parameter:
1. **Pre-fork (PHASE=1)**: Setup (zero DP, aliases, CC_VELOCITY_FLUX if CC_IBM)
2. **Branch B kernel (PHASE=2)**: PNV + Phases 2A-4b (species diffusion, heat diffusion, thermal) — called after WALL_BC in WORK_BRANCH=2
3. **Post-join kernel (PHASE=3)**: Phase 5+ (enthalpy advection, RTRM, species advection, sources, zone sums) — runs in WORK_BRANCH=2, copies WORK1_B→WORK1 for DIV_P2

**PNV ordering**: PREDICT_NORMAL_VELOCITY must run AFTER WALL_BC (WALL_BC sets B1%U_NORMAL_S for HVAC/SPECIFIED_MASS_FLUX walls). Placed at start of PHASE=2 so pipelined path (WALL_BC→PHASE=2) preserves correct ordering.

Requires Fortran kernel changes (splitting DIV_P1 into sections), unlike the simple VFLUX || WALL_BC fork.

**WORK conflict**: VFLUX uses WORK1-6, DIV_P1 early phases use WORK1-7. Requires per-branch scratch duplication (reuses Phase 2 infrastructure). WALL_BC uses no WORK arrays.

### DIV_P1 Phase Split (Predictor Only)

| Section | Phases | Ops/cell | UU/VV/WW? | FVX? | WALL_BC? |
|---------|--------|----------|-----------|------|----------|
| Pre-fork | Setup (zero DP, aliases, CC_VELOCITY_FLUX if CC_IBM) | ~2 | No | No | No |
| Branch B | PNV + 2A+2B (diffusion) + 3+4 (thermal) | 300 | **No** | **No** | Yes (PNV after WALL_BC) |
| Post-join | 5A-G (advection, sources) + 6 (zone sums) | ~1000 | **Read** | **Read** (CC_IBM) | Yes (wall corr.) |

### Checklist (Implementation)

- [x] Investigate Option B feasibility — full interior/wall split (result: NOT feasible)
- [x] Investigate DIV_P1 phase-level UU/VV/WW dependency (result: Phases 2A-4b are safe)
- [x] Verify CC_IBM safety for early phases (result: no special cases needed)
- [x] Split DIV_P1 kernel into pre-fork / early / late via PHASE parameter
- [x] Create C wrappers for split kernels (prefork, early_b, late_b)
- [x] Create predictor fork/join data types (pred_fork_data.h)
- [x] Create predictor fork state (PredForkState: MeshData → PredForkVFluxWork + PredForkDivWork)
- [x] Create predictor join state (PredJoinState: per-mesh matching)
- [x] Adapt VelocityFlux sub-graph for PredForkVFluxWork input (unwrap/wrap in pred_fork_vflux_subgraph.h)
- [x] Chain PARTICLE_MOMENTUM after VFLUX in Branch A
- [x] Chain DIV_P1_early after WALL_BC in Branch B (pred_fork_div_subgraph.h)
- [x] Wire predictor sub-graph with fork/join (CC_IBM sequential fallback)
- [x] Run verification suite (46/58 pass — no regressions from predictor pipelining)
- [x] Compare performance: VFLUX+PMOM hidden behind WBC+DIV — dancing_eddies: 1.12x (80ms hidden), fire_const_gamma: 1.10x (284ms hidden), species_props: 1.75x

---

## Performance Results

Measured with dot file execution stats on 3 test cases:

| Fork | dancing_eddies_4mesh | fire_const_gamma_2mesh | species_props_5mesh |
|------|---------------------|----------------------|-------------------|
| Phase 4: VFLUX \|\| COMB | 1.01x (no combustion) | **1.25x** (447ms hidden) | 1.11x |
| Phase 5: RAD \|\| DIV_P1 | **1.91x** (506ms hidden) | **1.29x** (534ms hidden) | 1.11x |
| Phase 6: VFLUX+PMOM \|\| WBC+DIV | 1.12x (80ms hidden) | 1.10x (284ms hidden) | **1.75x** |
| Combined savings/timestep | 587ms | 1264ms | 3.9ms |

**Key findings**:
- Fork 2 (RAD \|\| DIV_P1) is the highest-impact change — near-balanced branches give up to 1.91x on non-fire cases
- Fork 1 (VFLUX \|\| COMB) only helps when combustion is active — fire cases see 1.25x
- Predictor fork hides VFLUX+PMOM behind the larger WBC+DIV_P1_early critical path
- Branch imbalance limits speedup: WallBCBlockOrch dominates predictor Branch B (sequential preprocessing)
- All fork overheads (dispatch + join) are negligible (<10ms)

---

## Phase 7: Sequential Node Merging

**Objective**: Reduce Hedgehog overhead by merging sequential collector+barrier pairs into single BarrierState nodes. Add routine profiling and listing via `extraPrintingInformation()`.

### Changes

**New file**: `state/barrier_state.h` — Generic `BarrierState` + `BarrierStateManager` template
- Collects N MeshData, runs a `std::function` callback, re-emits N MeshData
- Built-in timing accumulator (total + average per invocation)
- `extraPrintingInformation()` shows routine names and timing stats in dot files

**Predictor sub-graph** — 5 collector+barrier pairs merged:
- [x] Collector(1) + MeshExchange(1) → BarrierState
- [x] PredHvacCollector + HvacInitDivTask → BarrierState
- [x] PredDivCollector + DivergenceExchangeTask → BarrierState
- [x] PredPressureCollector + PressureIterationTask → BarrierState (sequential pressure path)
- [x] Collector(3) + MeshExchange(3) → BarrierState

**Corrector sub-graph** — 7 collector+barrier pairs merged:
- [x] Collector(4) + MeshExchange(4) → BarrierState
- [x] SootHvacCollector + SootHvacTask → BarrierState
- [x] ParticleRemoveMoveCollector + RemoveMoveParticlesTask → BarrierState
- [x] Collector(7) + MeshExchange(7) → BarrierState
- [x] Collector(6a) + MeshExchange(6a) → BarrierState
- [x] CorrDivCollector + DivergenceExchangeTask → BarrierState
- [x] CorrPressureCollector + PressureIterationTask → BarrierState (sequential pressure path)
- [x] Collector(6b) + MeshExchange(6b) → BarrierState

**Main graph** — 3 nodes merged into 1:
- [x] TimestepGlobalTask + DumpMeshOutputsTask + TimestepDumpCollector → TimestepDumpState

**Kept as-is** (no collector to merge):
- PhaseTransitionTask (receives BarrierData from PredFinal sub-graph)
- MeshExchange(2) (receives BarrierData from CorrRadiation/Join2)
- ChangeTimeStepCollector (outputs BarrierData for ChangeTimeStep sub-graph input)
- PredPressureCollector / CorrPressureCollector (parallel pressure sub-graph input)

**Cleanup**:
- [x] Removed 10 barrier task classes from `barrier_tasks.h`
- [x] Added `extraPrintingInformation()` to remaining PhaseTransitionTask and MeshExchangeTask
- [x] `timestep_tasks.h` now dead code (TimestepGlobalTask + DumpMeshOutputsTask absorbed)

### Results
- **Node reduction**: ~16 fewer graph nodes, ~16 fewer inter-node queues
- **Custom tests**: 12/12 pass
- **Verification suite**: 46/58 pass (no regressions)
- **Dot files**: All merged nodes show routine names + timing stats

### Performance Analysis (dot file stats)

Parser: `test_cases/parse_dot_stats.py` — extracts per-node timing from Hedgehog dot files.

**Cross-test overhead comparison** (from `parse_dot_stats.py` on 16 test runs):

| Test Case | Graph(s) | Barrier(ms) | Orch(ms) | Coll(ms) | Overhead% | Kernel% |
|---|---|---|---|---|---|---|
| activate_sprinklers (1M) | 22.8 | 2064 | 1196 | 907 | 18.3% | 81.7% |
| bucket_test_1_short (1M) | 21.3 | 1214 | 1766 | 1366 | 22.9% | 77.1% |
| dancing_eddies_4mesh | 5.3 | 598 | 545 | 288 | 26.8% | 73.2% |
| dancing_eddies_ulmat | 16.4 | 2258 | 1636 | 885 | 29.2% | 70.8% |
| fire_const_gamma_2mesh | 30.9 | 1660 | 1406 | 3017 | 24.2% | 75.8% |
| species_props_5mesh | 0.1 | 14 | 2 | 3 | 18.6% | 81.4% |

**Key findings**:
- Barrier states account for 5-14% of graph time (dominated by Pressure + TimestepDump + RemoveMove)
- Orchestrator dispatch overhead: 4-10% (dominated by WallBCBlockOrch sequential work)
- Collector gather overhead: 4-10% (dominated by CorrFinalBlockColl and DensityCollector)
- Total state overhead: 18-29% → 71-82% of time in parallel kernel execution
- Merging collector+barrier eliminated ~16 inter-node queues per graph

**Block vs NoBlock vs Sequential** (bucket_test_1_short, single mesh):

| Variant | Graph(s) | Overhead% | Speedup |
|---|---|---|---|
| block (merged barriers) | 21.3 | 22.9% | 1.54x |
| noblock | 24.3 | 20.2% | 1.35x |
| sequential (1 thread) | 32.9 | 14.4% | 1.00x |

**Top overhead hotspots** (activate_sprinklers, 1 mesh, 1807 timesteps):
1. RemoveMove (barrier): 886ms — REMOVE_PARTICLES + MOVE_PARTICLES
2. WallBCBlockOrch (orchestrator): 501ms — sequential wall setup before parallel kernels
3. TimestepDump (barrier): 395ms — all dump/diagnostic I/O
4. CorrPressure (barrier): 377ms — pressure iteration
5. CorrFinalBlockColl (collector): 373ms — corrector final gather

---

## Summary

| Phase | Description | Fortran Changes | Graph Changes | WORK Branch | Expected Speedup |
|-------|-------------|-----------------|---------------|-------------|-----------------|
| 1 | Kernel extraction | SKIP_QR flag + QR kernel | None | No | N/A (validation) |
| 2 | Scratch arrays | WORK_B arrays + WORK_BRANCH | None | Infrastructure | N/A (validation) |
| 3 | Sequential driver | None | Task internals only | Validation | N/A (validation) |
| 4 | Corrector Fork 1 | None | Fork/join states | Not needed | ~2% corrector |
| 5 | Corrector Fork 2 | From Phase 1 | Fork/join states | Active | ~30% corrector |
| 6 | Predictor pipeline | DIV_P1 split into early/late | Fork/join states | Not needed | ~18% predictor |

**Phases 1-3**: Foundation work. No speedup, but validates all kernels.
**Phase 4**: First graph restructuring. Simple, low risk. Small but reliable speedup.
**Phase 5**: High-impact change. Relies on Phases 1-3 infrastructure.
**Phase 6**: Predictor fork with DIV_P1 early phases in Branch B. 266 ops/cell saved (4× improvement over simple VFLUX || WALL_BC).

---

## Test Criteria

Each phase must pass:
1. **Build**: `cmake --build . --target fds_hh -j$(nproc)` succeeds
2. **Custom tests**: `cd test_cases && python3 run_tests.py -v` (12 cases pass)
3. **Verification suite**: `cd test_cases && python3 run_verification.py test --no-redundant --max-gold-time 30 --timeout 120 --tolerance 1e-6` (46+ cases pass at 1e-6 tolerance)
4. **Dot file**: `graph->createDotFile(...)` shows correct pipeline structure
