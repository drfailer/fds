# WALL_BC Thread-Safe Conversions - Summary

## Successfully Completed (2 of 3 Quick Wins)

### ✅ 1. CALC_HVAC_BC (wall.f90)

**Approach**: Added explicit mesh and predictor flag arguments.

**Changes**:
```fortran
! Old signature:
SUBROUTINE CALC_HVAC_BC(BC, B1, SF)

! New signature:
SUBROUTINE CALC_HVAC_BC(M, PREDICTOR_FLAG, BC, B1, SF)
  TYPE(MESH_TYPE), INTENT(INOUT) :: M
  LOGICAL, INTENT(IN) :: PREDICTOR_FLAG
```

**Pointer replacement**:
```fortran
! Old (module-level pointer):
B1%RHO_F = PBAR_P(BC%KK, B1%PRESSURE_ZONE) / (RSUM_F*B1%TMP_G)

! New (conditional mesh access):
IF (PREDICTOR_FLAG) THEN
   B1%RHO_F = M%PBAR_S(BC%KK, B1%PRESSURE_ZONE) / (RSUM_F*B1%TMP_G)
ELSE
   B1%RHO_F = M%PBAR(BC%KK, B1%PRESSURE_ZONE) / (RSUM_F*B1%TMP_G)
ENDIF
```

**Call sites updated**: 2 locations in WALL_BC

---

### ✅ 2. HEAT_TRANSFER_COEFFICIENT (func.f90)

**Approach**: Use integer indices instead of pointers to avoid Fortran ALLOCATABLE/TARGET limitation.

**Key insight**: Instead of creating pointer variables to array elements, store the indices and use direct array access.

**Changes**:
```fortran
! Old signature:
REAL(EB) FUNCTION HEAT_TRANSFER_COEFFICIENT(NM, T, ...)
  INTEGER, INTENT(IN) :: NM
  TYPE(MESH_TYPE), POINTER :: M
  TYPE(WALL_TYPE), POINTER :: WC
  TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
  ...
  M => MESHES(NM)
  WC => M%WALL(WALL_INDEX_IN)
  B1 => M%BOUNDARY_PROP1(WC%B1_INDEX)

! New signature:
REAL(EB) FUNCTION HEAT_TRANSFER_COEFFICIENT(M, T, ...)
  TYPE(MESH_TYPE), INTENT(INOUT) :: M
  INTEGER :: B1_INDEX, B2_INDEX, BC_INDEX
  ...
  B1_INDEX = M%WALL(WALL_INDEX_IN)%B1_INDEX
  B2_INDEX = M%WALL(WALL_INDEX_IN)%B2_INDEX
  BC_INDEX = M%WALL(WALL_INDEX_IN)%BC_INDEX
```

**Access pattern replacement** (automated with awk):
```fortran
! Old:
B1%TMP_F

! New:
M%BOUNDARY_PROP1(B1_INDEX)%TMP_F
```

**Additional fixes**:
- Broke 4 lines exceeding 132-character Fortran limit
- Updated 15+ call sites in wall.f90 and dump.f90

---

### ⏸️ 3. DEPOSIT_PARTICLE_MASS (wall.f90)

**Status**: Already thread-safe, no conversion needed.

**Rationale**: Does not use module-level pointers (PBAR_P, RHOP, UU, VV, WW, ZZP).

**Future**: Could be moved to wall_kernels.f90 using the same index-based approach as HEAT_TRANSFER_COEFFICIENT, but it's not urgent since it's already safe in its current location.

---

## Key Technical Solution: Index-Based Access

**Problem**: Fortran gfortran doesn't allow pointer assignment to elements of ALLOCATABLE arrays unless they have the TARGET attribute:
```fortran
TYPE(MESH_TYPE), INTENT(INOUT) :: M
TYPE(WALL_TYPE), POINTER :: WC
WC => M%WALL(WALL_INDEX)  ! Error: target is neither TARGET nor POINTER
```

**Solution**: Don't use pointers—use integer indices and direct array access:
```fortran
TYPE(MESH_TYPE), INTENT(INOUT) :: M
INTEGER :: WC_B1_INDEX
WC_B1_INDEX = M%WALL(WALL_INDEX)%B1_INDEX
...
M%BOUNDARY_PROP1(WC_B1_INDEX)%HEAT_TRANS_COEF = ...
```

**Trade-off**: More verbose code, but:
- ✅ Thread-safe (no module-level state)
- ✅ No need to modify MESH_TYPE definition
- ✅ Works with existing Fortran compiler
- ✅ Explicit mesh parameter enables parallelization

---

## Build Verification

✅ **All changes compile successfully**
- Target: `fds`
- Build system: CMake
- Compiler: gfortran
- Warnings: None (all treated as errors)

---

## Summary Table

| Routine | Status | Thread-Safe? | Approach | Call Sites Updated |
|---------|--------|--------------|----------|-------------------|
| CALC_HVAC_BC | ✅ Converted | ✅ Yes | Explicit M + PREDICTOR_FLAG | 2 |
| HEAT_TRANSFER_COEFFICIENT | ✅ Converted | ✅ Yes | Index-based access (no pointers) | 15+ |
| DEPOSIT_PARTICLE_MASS | ⏸️ Unchanged | ✅ Yes | Already safe (no module pointers) | 0 |

**Net result**: 2 of 3 converted, all 3 are thread-safe.

---

## Files Modified

1. **Source/func.f90**:
   - HEAT_TRANSFER_COEFFICIENT signature and implementation
   - Replaced all pointer variables with index-based access
   - Added line continuations for long lines

2. **Source/wall.f90**:
   - CALC_HVAC_BC signature and implementation
   - Updated all HEAT_TRANSFER_COEFFICIENT call sites
   - Updated all CALC_HVAC_BC call sites
   - Added line continuation for long HEAT_TRANSFER_COEFFICIENT call

3. **Source/dump.f90**:
   - Updated HEAT_TRANSFER_COEFFICIENT call sites

---

## Next Steps for Full WALL_BC Parallelization

With the quick wins complete, the critical path is now:

### Phase A (Estimated 2-3 hours each):
1. **SURFACE_HEAT_TRANSFER** - Largest impact, 90% of wall cells
   - 379 lines, uses UU, VV, WW, RHOP, ZZP, PBAR_P extensively
   - Split INTERPOLATED_BC case to separate routine
   - Convert to index-based access (same pattern as HEAT_TRANSFER_COEFFICIENT)

2. **CALCULATE_ZZ_F** - Called for all wall cells
   - 414 lines, uses UU, VV, WW for normal velocity
   - Extract CONSUME_MASS section to separate routine
   - Convert remaining to index-based access

### Phase B (Lower priority, 3-4 hours):
3. **SOLID_HEAT_TRANSFER** - Much work stays in Phase 3 anyway
   - 1337 lines, complex heat conduction solver
   - BACK_MESH cases must remain sequential
   - Lower urgency since many cells deferred to Phase 3

**Total estimated effort to unblock parallelization**: 4-6 hours (Phases A only)

---

## Lessons Learned

1. **Avoid pointers to ALLOCATABLE elements**: Use indices instead
2. **Automated replacements save time**: awk script replaced 17 `B1%` → `M%BOUNDARY_PROP1(B1_INDEX)%`
3. **Watch line length**: Fortran 132-char limit catches verbose array access
4. **Test incrementally**: Build after each conversion to catch issues early
5. **User feedback is valuable**: "Why use WC pointer?" led to simpler solution

This approach (index-based access) can be applied to all remaining routines without modifying MESH_TYPE or adding TARGET attributes.
