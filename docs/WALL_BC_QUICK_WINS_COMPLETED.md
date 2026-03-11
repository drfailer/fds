# WALL_BC "Quick Wins" - Completion Summary

## Completed Conversions

### ✅ CALC_HVAC_BC - Successfully Converted

**File**: `Source/wall.f90`

**Changes**:
1. **Signature updated** to accept explicit mesh and predictor flag:
   ```fortran
   SUBROUTINE CALC_HVAC_BC(M, PREDICTOR_FLAG, BC, B1, SF)
     TYPE(MESH_TYPE), INTENT(INOUT) :: M
     LOGICAL, INTENT(IN) :: PREDICTOR_FLAG
   ```

2. **Module-level pointer replaced** with conditional mesh access:
   ```fortran
   ! Old:
   B1%RHO_F = PBAR_P(BC%KK,B1%PRESSURE_ZONE) / (RSUM_F*B1%TMP_G)

   ! New:
   IF (PREDICTOR_FLAG) THEN
      B1%RHO_F = M%PBAR_S(BC%KK,B1%PRESSURE_ZONE) / (RSUM_F*B1%TMP_G)
   ELSE
      B1%RHO_F = M%PBAR(BC%KK,B1%PRESSURE_ZONE) / (RSUM_F*B1%TMP_G)
   ENDIF
   ```

3. **Call sites updated** (2 locations in WALL_BC):
   ```fortran
   ! Old:
   CALL CALC_HVAC_BC(BC, B1, SF)

   ! New:
   CALL CALC_HVAC_BC(MESHES(NM), PREDICTOR, BC, B1, SF)
   ```

**Status**: ✅ Compiles successfully, thread-safe

---

## Deferred Conversions

### ⏸️ HEAT_TRANSFER_COEFFICIENT - Deferred

**File**: `Source/func.f90`

**Reason for deferral**: Fortran compiler limitation with pointer assignment to ALLOCATABLE array components.

**Issue**: When `M` is passed as `INTENT(INOUT)`, the compiler doesn't allow pointer assignments like:
```fortran
WC => M%WALL(WALL_INDEX_IN)      ! Error: target is neither TARGET nor POINTER
B1 => M%BOUNDARY_PROP1(WC%B1_INDEX)
```

**Workaround options**:
1. Add `TARGET` attribute to ALLOCATABLE arrays in MESH_TYPE (requires changes to mesh.f90)
2. Refactor to use direct array access instead of pointers (extensive changes ~300 lines)
3. Pass B1, B2, BC as arguments instead of looking them up inside the function (signature change)

**Impact**: HEAT_TRANSFER_COEFFICIENT is called by SURFACE_HEAT_TRANSFER, which is a larger conversion target anyway. We'll tackle both together in the next phase.

---

### ⏸️ DEPOSIT_PARTICLE_MASS - Already Thread-Safe, Not Moved

**File**: `Source/wall.f90` (remains here)

**Status**: Already thread-safe (doesn't use module-level pointers)

**Reason for not moving**: Encountered same pointer assignment issues as HEAT_TRANSFER_COEFFICIENT. Since it's already thread-safe, leaving it in wall.f90 is acceptable.

**Future**: Can be moved to wall_kernels.f90 later if we add TARGET attributes to MESH_TYPE arrays.

---

## Build Verification

✅ **Build successful**: `cmake --build . --target fds` completes without errors

---

## Summary

| Routine | Status | Thread-Safe? | In Kernel Module? |
|---------|--------|--------------|-------------------|
| CALC_HVAC_BC | ✅ Converted | ✅ Yes | No (in wall.f90) |
| HEAT_TRANSFER_COEFFICIENT | ⏸️ Deferred | ❌ Uses MESHES(NM) | No (in func.f90) |
| DEPOSIT_PARTICLE_MASS | ⏸️ Not moved | ✅ Yes | No (in wall.f90) |

**Net result**: 1 of 3 routines converted to thread-safe kernel pattern.

---

## Next Steps

To complete the "quick wins" and enable further progress, we need to tackle the pointer assignment issue. Recommended approach:

**Option A (Cleanest)**: Add `TARGET` attribute to allocatable arrays in MESH_TYPE
```fortran
! In mesh.f90, MESH_TYPE definition:
TYPE(WALL_TYPE), ALLOCATABLE, TARGET, DIMENSION(:) :: WALL
TYPE(BOUNDARY_PROP1_TYPE), ALLOCATABLE, TARGET, DIMENSION(:) :: BOUNDARY_PROP1
TYPE(BOUNDARY_COORD_TYPE), ALLOCATABLE, TARGET, DIMENSION(:) :: BOUNDARY_COORD
...
```

This would allow pointer assignment to work and unblock:
- HEAT_TRANSFER_COEFFICIENT conversion
- DEPOSIT_PARTICLE_MASS move to wall_kernels.f90
- Future SURFACE_HEAT_TRANSFER conversion

**Estimated effort**: 30 minutes (modify MESH_TYPE, rebuild, test)

---

## Files Modified

1. `Source/wall.f90`:
   - CALC_HVAC_BC signature and implementation updated
   - Call sites updated to pass MESHES(NM) and PREDICTOR

2. `Source/type.f90`:
   - Added HAS_INTERPOLATED_BC and HAS_BACK_MESH flags to WALL_TYPE

3. `Source/init.f90`:
   - Added flag initialization in FIND_WALL_BACK_INDICES

**Build status**: ✅ All changes compile successfully
