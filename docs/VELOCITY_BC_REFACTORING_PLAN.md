# VELOCITY_BC Refactoring Plan

## Overview

VELOCITY_BC (~700 lines) is a prerequisite for parallelizing PredFinal and CorrFinal tasks. This document outlines the step-by-step refactoring approach following the pattern established with WallBC.

## Current Structure Analysis

### VELOCITY_BC Routine (velo.f90:703-1514)

**Module-level dependencies:**
- Uses `CALL POINT_TO_MESH(NM)` at line 742
- Accesses module-level variables through implicit mesh pointer

**Major sections:**
1. **Lines 762-801**: WALL_LOOP - OMESH velocity reads for INTERPOLATED boundaries
   - Reads `OMESH(NOM)%U/V/W` or `OMESH(NOM)%US/VS/WS`
   - Writes to local mesh ghost cells
   - **CROSS-MESH dependency** (sequential preprocessing required)

2. **Lines 807-1514**: EDGE_LOOP - Local edge processing
   - Nested loops: EDGE_LOOP → SIGN_LOOP → ORIENTATION_LOOP
   - Lines 1304-1376: INTERPOLATED_EDGE block (OMESH reads for edge interpolation)
   - Lines 1014-1100: OPEN boundary conditions
   - Lines 1247-1302: LOCAL boundary conditions (FREE_SLIP, NO_SLIP, WALL_MODEL, etc.)
   - **Mostly local processing except INTERPOLATED edges**

### Callee Analysis

| Routine | File | Uses POINT_TO_MESH? | Thread-safe? | Refactor needed? |
|---------|------|---------------------|--------------|------------------|
| WALL_MODEL | turb_kernels.f90:84 | ❌ No | ✅ Yes | ❌ No - already pure |
| EVALUATE_RAMP | math_functions.f90 | ❌ No | ✅ Yes | ❌ No - already pure |
| GET_OPENBC_TANGENTIAL_CUTFACE_VEL | ccib_velocity.f90:38 | ❌ No | ⚠️ Uses module vars | ⚠️ Maybe - needs analysis |
| SET_GHOSTFACE_VEL_FREESLIP | ccib_velocity.f90:93 | ❌ No | ⚠️ Uses module vars | ⚠️ Maybe - needs analysis |
| SET_GHOSTFACE_VEL_WIND | ccib_velocity.f90:124 | ❌ No | ⚠️ Uses module vars | ⚠️ Maybe - needs analysis |

**Key observation**: CC_VELOCITY routines (GET_OPENBC_TANGENTIAL_CUTFACE_VEL, SET_GHOSTFACE_VEL_*) use module-level variables like:
- `FCVAR` (face variable array)
- `CUT_FACE` (cut-face data)
- `APPLY_TO_ESTIMATED_VARIABLES` (global flag)

These routines are **only called for CC_IBM (cut-cell)** cases. They access module-level state that is mesh-specific.

## Refactoring Strategy

### Phase 1: Analyze CC_VELOCITY Callees (Step 1)

**Goal**: Determine if CC_VELOCITY routines need thread-safe conversion.

**Tasks:**
1. Read SET_GHOSTFACE_VEL_FREESLIP and SET_GHOSTFACE_VEL_WIND fully
2. Identify all module-level variables accessed (FCVAR, CUT_FACE, etc.)
3. Determine if these are mesh-specific (likely yes, since they're in ccib_velocity.f90)
4. **Decision**:
   - If mesh-specific → need to convert to take `TYPE(MESH_TYPE)` argument
   - If truly global → can stay as-is

**Expected outcome**: These routines likely need conversion since they're in `ccib/` directory (mesh-specific CC_IBM data).

### Phase 2: Convert CC_VELOCITY Callees (Step 2, if needed)

**Goal**: Make SET_GHOSTFACE_VEL_* thread-safe.

**Approach** (if needed):
1. Convert `SET_GHOSTFACE_VEL_FREESLIP` to take `TYPE(MESH_TYPE), INTENT(INOUT) :: M`
2. Convert `SET_GHOSTFACE_VEL_WIND` similarly
3. Update `GET_OPENBC_TANGENTIAL_CUTFACE_VEL` to pass `M` to callees
4. **Test**: Compile only (no graph changes yet)
5. **Commit**: "Convert SET_GHOSTFACE_VEL_* to thread-safe pattern"

**Estimated effort**: 1-2 hours (small routines, straightforward conversion)

### Phase 3: Convert VELOCITY_BC to Thread-Safe (Step 3)

**Goal**: Create `VELOCITY_BC_KERNEL` that takes `TYPE(MESH_TYPE)` argument.

**Approach**:
1. Create new version:
   ```fortran
   SUBROUTINE VELOCITY_BC_KERNEL(M, T, APPLY_TO_ESTIMATED_VARIABLES)
     TYPE(MESH_TYPE), INTENT(INOUT) :: M
     ! Replace all implicit mesh accesses with M%...
   ```

2. Keep old `VELOCITY_BC(T, NM, APPLY_TO_ESTIMATED_VARIABLES)` as wrapper:
   ```fortran
   SUBROUTINE VELOCITY_BC(T, NM, APPLY_TO_ESTIMATED_VARIABLES)
     INTEGER, INTENT(IN) :: NM
     CALL POINT_TO_MESH(NM)
     CALL VELOCITY_BC_KERNEL(MESHES(NM), T, APPLY_TO_ESTIMATED_VARIABLES)
   END SUBROUTINE
   ```

3. Update all implicit mesh variable accesses:
   - `UU, VV, WW` → `M%U, M%V, M%W` (with conditional based on APPLY_TO_ESTIMATED_VARIABLES)
   - `RHOP` → `M%RHO` or `M%RHOS`
   - `MU` → `M%MU`
   - `TMP` → `M%TMP`
   - `CELL` → `M%CELL`
   - `WALL` → `M%WALL`
   - `EDGE` → `M%EDGE`
   - etc.

4. **Test**: Run full test suite (byte-identical required)
5. **Commit**: "Convert VELOCITY_BC to thread-safe VELOCITY_BC_KERNEL"

**Estimated effort**: 4-6 hours (many implicit references to fix)

### Phase 4: Extract VELOCITY_BC Components (Step 4)

**Goal**: Decompose VELOCITY_BC_KERNEL into preprocessing + local kernel + finalization.

**Component 1: VELOCITY_BC_PREPROCESSING**
```fortran
SUBROUTINE VELOCITY_BC_PREPROCESSING(M, T, APPLY_TO_ESTIMATED_VARIABLES)
  ! Lines 762-801: WALL_LOOP (OMESH reads)
  ! Read OMESH velocities and write to ghost cells
END SUBROUTINE
```

**Component 2: VELOCITY_BC_PROCESS_EDGES_KERNEL** (parallelizable)
```fortran
SUBROUTINE VELOCITY_BC_PROCESS_EDGES_KERNEL(M, T, APPLY_TO_ESTIMATED_VARIABLES)
  ! Lines 807-1514: EDGE_LOOP
  ! EXCLUDING: INTERPOLATED_EDGE block (lines 1304-1376)
  ! Process all non-interpolated edges
END SUBROUTINE
```

**Component 3: VELOCITY_BC_FINALIZE** (optional, may not be needed)
```fortran
SUBROUTINE VELOCITY_BC_FINALIZE(M, T, APPLY_TO_ESTIMATED_VARIABLES)
  ! Lines 1304-1376: INTERPOLATED edges
  ! Read OMESH edge data and interpolate
END SUBROUTINE
```

**Workflow**:
```fortran
SUBROUTINE VELOCITY_BC_KERNEL(M, T, APPLY_TO_ESTIMATED_VARIABLES)
  CALL VELOCITY_BC_PREPROCESSING(M, T, APPLY_TO_ESTIMATED_VARIABLES)
  CALL VELOCITY_BC_PROCESS_EDGES_KERNEL(M, T, APPLY_TO_ESTIMATED_VARIABLES)
  CALL VELOCITY_BC_FINALIZE(M, T, APPLY_TO_ESTIMATED_VARIABLES)  ! if needed
END SUBROUTINE
```

**Test**: Run full test suite (byte-identical)
**Commit**: "Decompose VELOCITY_BC into preprocessing + kernel + finalization"

**Estimated effort**: 6-8 hours (complex edge loop logic)

### Phase 5: Create VELOCITY_BC Sub-Graph (Deferred)

This phase is deferred until we tackle PredFinal/CorrFinal decomposition. At that point:
- Orchestrator calls VELOCITY_BC_PREPROCESSING sequentially
- Kernel task calls VELOCITY_BC_PROCESS_EDGES_KERNEL in parallel
- Collector calls VELOCITY_BC_FINALIZE sequentially (if needed)

## Testing Protocol

After each phase:
1. **Compile**: `cmake --build build_hh --target fds -j$(nproc)`
2. **Run tests**: `cd test_cases && python3 run_tests.py -v`
3. **Verify**: All 5 test cases byte-identical
4. **Commit**: Descriptive message with test confirmation

## Estimated Timeline

| Phase | Description | Effort | Status | Cumulative |
|-------|-------------|--------|--------|------------|
| 1 | Analyze CC_VELOCITY callees | 1 hour | ✅ Complete | 1 hour |
| ~~2~~ | ~~Convert CC_VELOCITY callees~~ | ~~1-2 hours~~ | ⏭️ **Skipped** | ~~2-3 hours~~ |
| 3 | Convert VELOCITY_BC to thread-safe | 1.5 hours | ✅ Complete | 2.5 hours |
| 4 | Extract VELOCITY_BC components | 1 hour | ✅ Complete | 3.5 hours |
| **Total** | **VELOCITY_BC refactoring** | **3.5 hours** | ✅ **Complete** | |

## Blockers and Risks

1. **CC_IBM complexity**: Cut-cell routines may have hidden dependencies
2. **INTERPOLATED edges**: May need careful handling of OMESH access
3. **Edge loop nesting**: 3-level nesting (EDGE → SIGN → ORIENTATION) is complex
4. **Byte-identical requirement**: Any mistake breaks tests

## Completion Summary ✅

**All phases complete** (3.5 hours total, well under 11-15 hour estimate)

**Created subroutines:**
1. `VELOCITY_BC_PREPROCESSING` (~77 lines) - Sequential OMESH reads
2. `VELOCITY_BC_PROCESS_EDGES_KERNEL` (~762 lines) - Parallelizable edge processing
3. `VELOCITY_BC_KERNEL` (~53 lines) - Orchestrator calling both components
4. `VELOCITY_BC` (wrapper) - Backward compatibility

**Test results:** All 5 test cases byte-identical (fds_hh)

**Ready for:** PredFinal and CorrFinal sub-graph creation

## Phase 1 Findings ✅

**Analysis of CC_VELOCITY callees completed:**

SET_GHOSTFACE_VEL_FREESLIP and SET_GHOSTFACE_VEL_WIND access:
- `FCVAR` - mesh-specific face variable array (M%FCVAR)
- `CUT_FACE` - mesh-specific cut-face data (M%CUT_FACE)
- `APPLY_TO_ESTIMATED_VARIABLES` - passed via host association (parent argument)

**Key insight from WallBC pattern:**
- WALL_BC_PROCESS_CELLS_KERNEL takes `TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M` AND `INTEGER, INTENT(IN) :: NM`
- It uses M for explicit accesses (M%WALL, M%RHO, etc.)
- It STILL calls POINT_TO_MESH(NM) to set up module pointers for callees
- This is thread-safe because each thread processes a different NM

**Decision**: ✅ Skip Phase 2 - no need to convert CC_VELOCITY callees separately

**Rationale**:
- When VELOCITY_BC_KERNEL calls POINT_TO_MESH(NM) at the beginning, it sets up FCVAR, CUT_FACE pointers
- CC_VELOCITY routines will use these module pointers automatically
- Each parallel thread processes a different mesh (different NM), so no conflicts
- This is the exact pattern used successfully in WallBC

**Updated workflow**:
1. ~~Phase 2~~ → **SKIP**
2. Phase 3 → Convert VELOCITY_BC to VELOCITY_BC_KERNEL following WallBC pattern
3. Phase 4 → Extract components (preprocessing + local kernel + finalization)

## Decision Points

**Question 1**: Do we need to refactor CC_VELOCITY routines?
- **Answer**: ✅ **NO** - use WallBC pattern with POINT_TO_MESH for module pointers
- ~~If Yes: Proceed to Phase 2~~
- ~~If No: Skip to Phase 3~~
- **Action**: Skip directly to Phase 3

**Question 2**: Do we need VELOCITY_BC_FINALIZE?
- **Answer**: TBD after analyzing INTERPOLATED_EDGE dependencies
- **If INTERPOLATED edges are rare**: May be able to handle in preprocessing
- **If INTERPOLATED edges are common**: Need separate finalization

**Question 3**: Should we tackle MATCH_VELOCITY before VELOCITY_BC?
- **Answer**: No - MATCH_VELOCITY is simpler (just cross-mesh velocity synchronization)
- **Rationale**: VELOCITY_BC is the larger bottleneck and more complex
