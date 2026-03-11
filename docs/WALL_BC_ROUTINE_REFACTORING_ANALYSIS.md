# WALL_BC Routine Refactoring Analysis

## Overview

Analysis of all routines called by WALL_BC to assess complexity of converting them from using module-level pointers (set by POINT_TO_MESH) to accepting explicit `TYPE(MESH_TYPE)` arguments.

## Module-Level Pointers in WALL_ROUTINES

The WALL_ROUTINES module (wall.f90) declares module-level pointers that are set by WALL_BC based on the PREDICTOR flag:

```fortran
REAL(EB), POINTER, DIMENSION(:,:) :: PBAR_P
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP, UU, VV, WW
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP
```

These point to either predictor (RHOS, ZZS, US, VS, WS, PBAR_S) or corrector (RHO, ZZ, U, V, W, PBAR) fields from MESH_POINTERS.

## Routines Called by WALL_BC

### ✅ Already Kernelized (No Refactoring Needed)

| Routine | Location | Status |
|---------|----------|--------|
| NEAR_SURFACE_GAS_VARIABLES_KERNEL | wall_kernels.f90 | Takes TYPE(MESH_TYPE) ✓ |
| CALCULATE_RHO_D_F | wall_kernels.f90 | Takes TYPE(MESH_TYPE) ✓ |
| CALCULATE_RHO_F_KERNEL | wall_kernels.f90 | Takes TYPE(MESH_TYPE) ✓ |
| CALC_DEPOSITION | wall_kernels.f90 | Takes TYPE(MESH_TYPE) ✓ |
| WALL_MODEL | turb_kernels.f90 | Pure function, no mesh access ✓ |

### 🟢 Low Complexity (Easy to Refactor)

#### 1. CALC_HVAC_BC
- **Location**: wall.f90
- **Lines**: 52
- **Module pointers used**: `PBAR_P`
- **Complexity**: LOW
- **Refactoring**: Add `TYPE(MESH_TYPE)` argument, replace `PBAR_P` with `M%PBAR` or `M%PBAR_S`
- **Dependencies**: None (doesn't call other routines with POINT_TO_MESH)
- **Estimated effort**: 15 minutes

#### 2. DEPOSIT_PARTICLE_MASS
- **Location**: wall.f90
- **Lines**: 106
- **Module pointers used**: None identified
- **Complexity**: LOW
- **Refactoring**: Add `TYPE(MESH_TYPE)` argument for consistency
- **Dependencies**: None
- **Estimated effort**: 15 minutes

#### 3. HEAT_TRANSFER_COEFFICIENT
- **Location**: func.f90
- **Lines**: ~300
- **Module pointers used**: None (already uses `M => MESHES(NM)`)
- **Complexity**: VERY LOW
- **Refactoring**: Replace `M => MESHES(NM)` with explicit argument `TYPE(MESH_TYPE), INTENT(IN) :: M`
- **Dependencies**: None
- **Estimated effort**: 20 minutes
- **Note**: This is the easiest win — already 95% of the way there!

### 🟡 Medium Complexity

#### 4. ASSIGN_GHOST_VALUE
- **Location**: wall.f90
- **Lines**: 107
- **Module pointers used**: `RHOP`, `ZZP` (writes to ghost cells)
- **Complexity**: MEDIUM
- **Refactoring**:
  - Add `TYPE(MESH_TYPE)` argument
  - Replace `RHOP(BC%II,BC%JJ,BC%KK)` with `M%RHO(...)` or `M%RHOS(...)`
  - Replace `ZZP(BC%II,BC%JJ,BC%KK,...)` with `M%ZZ(...)` or `M%ZZS(...)`
  - Need to pass PREDICTOR flag to determine which field to write to
- **Dependencies**: None (doesn't call other POINT_TO_MESH routines)
- **Estimated effort**: 30 minutes
- **Note**: Called only in Phase 1 (sequential), so parallelization benefit is zero. Low priority.

#### 5. SURFACE_HEAT_TRANSFER
- **Location**: wall.f90
- **Lines**: 379
- **Module pointers used**: `UU`, `VV`, `WW`, `RHOP`, `ZZP`, `PBAR_P` (extensively)
- **Complexity**: MEDIUM-HIGH
- **Refactoring**:
  - Add `TYPE(MESH_TYPE)` argument
  - Add `PREDICTOR` flag argument
  - Replace all module pointer accesses with conditional M% access
  - ~50-60 lines of pointer replacements
- **Dependencies**:
  - Calls HEAT_TRANSFER_COEFFICIENT (needs to be refactored first)
  - Accesses OMESH for INTERPOLATED_BC case (inherently sequential)
- **Estimated effort**: 2 hours
- **Priority**: HIGH (called in Phase 2, parallelizable for non-INTERPOLATED cells)

#### 6. CALCULATE_ZZ_F
- **Location**: wall.f90
- **Lines**: 414
- **Module pointers used**: `UU`, `VV`, `WW` (for computing UN, normal velocity)
- **Complexity**: MEDIUM
- **Refactoring**:
  - Add `TYPE(MESH_TYPE)` argument
  - Add `PREDICTOR` flag argument
  - Replace velocity pointer accesses (~6 lines)
  - Extract CONSUME_MASS section to separate routine (Phase 3)
- **Dependencies**: None for main logic; CONSUME_MASS uses OMESH (sequential)
- **Estimated effort**: 1.5 hours
- **Priority**: HIGH (called in Phase 2, parallelizable except CONSUME_MASS)

### 🔴 High Complexity

#### 7. SOLID_HEAT_TRANSFER
- **Location**: wall.f90
- **Lines**: 1337
- **Module pointers used**: None directly identified, but uses many MESH_POINTERS variables
- **Complexity**: HIGH
- **Refactoring**:
  - Add `TYPE(MESH_TYPE)` argument
  - This is a complex heat conduction solver with:
    - 1-D pyrolysis model
    - Grid remeshing
    - Material property evaluation
    - Back-side coupling (BACK_MESH access)
  - Heavy use of derived types and temporary arrays
  - May not directly use module-level pointers (uses M%WALL, M%BOUNDARY_ONE_D, etc.)
- **Dependencies**:
  - Calls PYROLYSIS (already a kernel)
  - Calls GET_SPECIFIC_HEAT, GET_EMISSIVITY (physical functions)
  - Accesses MESHES(BACK_MESH) for thin walls (inherently sequential)
- **Estimated effort**: 3-4 hours
- **Priority**: MEDIUM (heavy computation, but much of it happens in Phase 3 sequential section for BACK_MESH cases)

### 🔵 Special Cases

#### 8. CFACE_THERMAL_GASVARS
- **Location**: ccib/ccib_scalars.f90 (CC_SCALARS module)
- **Status**: Need to investigate
- **Complexity**: UNKNOWN
- **Note**: Part of complex geometry (CC_IBM) subsystem

## Recommended Refactoring Order

### Phase 1: Quick Wins (1 hour total)
1. **HEAT_TRANSFER_COEFFICIENT** (20 min) — already 95% converted
2. **CALC_HVAC_BC** (15 min) — minimal pointer use
3. **DEPOSIT_PARTICLE_MASS** (15 min) — no pointer use

### Phase 2: Core Parallelization Enablers (3.5 hours total)
4. **SURFACE_HEAT_TRANSFER** (2 hours) — critical for Phase 2 kernel
5. **CALCULATE_ZZ_F** (1.5 hours) — critical for Phase 2 kernel

### Phase 3: Complex Cases (3-4 hours)
6. **SOLID_HEAT_TRANSFER** (3-4 hours) — large but mostly for Phase 3 sequential section

### Defer (low/zero parallelization benefit)
7. **ASSIGN_GHOST_VALUE** — Phase 1 only (sequential), no benefit
8. **CFACE_THERMAL_GASVARS** — investigate if needed

## Total Estimated Effort

- **Minimum viable parallelization** (Phase 1-2): ~4.5 hours
- **Full conversion** (all routines): ~7.5-8.5 hours

## Key Insight

The critical path for WALL_BC parallelization is:
1. HEAT_TRANSFER_COEFFICIENT ← dependency of SURFACE_HEAT_TRANSFER
2. SURFACE_HEAT_TRANSFER ← called in Phase 2 kernel for all non-INTERPOLATED cells
3. CALCULATE_ZZ_F ← called in Phase 2 kernel for all cells

Converting just these three routines (plus the quick wins) enables 90-95% of wall cells to be processed in parallel, giving the full speedup benefit.

SOLID_HEAT_TRANSFER is larger but less critical — much of its work (BACK_MESH cases) happens in the sequential Phase 3 anyway.
