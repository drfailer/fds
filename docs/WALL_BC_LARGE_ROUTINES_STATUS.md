# SURFACE_HEAT_TRANSFER and CALCULATE_ZZ_F Thread-Safe Conversion

## Status: In Progress

### Challenge: Fortran TARGET Attribute Limitation

**Problem**: Cannot add TARGET attribute to ALLOCATABLE arrays inside TYPE definitions.

```fortran
! This is NOT allowed in Fortran:
TYPE MESH_TYPE
   REAL(EB), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: U  ! ERROR
END TYPE
```

**Error**:
```
Error: Attribute at (1) is not allowed in a TYPE definition
```

**Impact**: Without TARGET, we cannot use pointer assignment to ALLOCATABLE array elements:
```fortran
TYPE(MESH_TYPE), INTENT(INOUT) :: M
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU
UU => M%U  ! ERROR: target is neither TARGET nor POINTER
```

### Required Approach: No-Pointer Direct Access

Since we cannot:
1. Add TARGET to MESH_TYPE arrays (language limitation)
2. Use pointers to ALLOCATABLE arrays without TARGET

We must convert to direct array access without any pointers.

### Conversion Pattern for Large Routines

**Before** (module-level pointers):
```fortran
! Module level (set by POINT_TO_MESH):
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU, VV, WW

! In routine:
UN = UU(II,JJ,KK)
```

**After** (conditional direct access):
```fortran
! In routine:
REAL(EB) :: UN
IF (PREDICTOR_FLAG) THEN
   UN = M%US(II,JJ,KK)
ELSE
   UN = M%U(II,JJ,KK)
ENDIF
```

**Verbosity trade-off**: Much more verbose but thread-safe.

### Routines Under Conversion

1. **SURFACE_HEAT_TRANSFER** (379 lines)
   - Status: Partially converted
   - Signature updated: `(M, PREDICTOR_FLAG, T, SF, BC, B1, ...)`
   - Array accesses updated: TMP, RSUM, MU, DX, DY, DZ, etc. → M%
   - Remaining: Remove pointer assignments, use direct conditional access

2. **CALCULATE_ZZ_F** (413 lines)
   - Status: Partially converted
   - Signature updated: `(M, PREDICTOR_FLAG, T, DT, ...)`
   - Array accesses updated: similar to SURFACE_HEAT_TRANSFER
   - Remaining: Remove pointer assignments, use direct conditional access

### Arrays Requiring Conditional Access

Based on PREDICTOR_FLAG:
- `RHOP`: `M%RHOS` (predictor) or `M%RHO` (corrector)
- `ZZP`: `M%ZZS` (predictor) or `M%ZZ` (corrector)
- `UU, VV, WW`: `M%US, M%VS, M%WS` (predictor) or `M%U, M%V, M%W` (corrector)
- `PBAR_P`: `M%PBAR_S` (predictor) or `M%PBAR` (corrector)

### Mesh Component Pointers

Cannot use:
```fortran
WC => M%WALL(WALL_INDEX)  ! ERROR: no TARGET
```

Must use:
```fortran
! Direct access:
IF (M%WALL(WALL_INDEX)%VENT_INDEX > 0) THEN
   ! Use M%VENTS(M%WALL(WALL_INDEX)%VENT_INDEX)%...
ENDIF
```

Or store indices:
```fortran
VENT_INDEX = M%WALL(WALL_INDEX)%VENT_INDEX
IF (VENT_INDEX > 0) THEN
   ! Use M%VENTS(VENT_INDEX)%...
ENDIF
```

### Estimated Effort

Given the size and complexity:
- ~100-150 pointer uses across both routines
- Each requires conditional or direct access replacement
- Estimated time: 4-6 hours for complete conversion
- Risk: High complexity, easy to introduce bugs

### Alternative: Smaller Decomposition

**Recommendation**: Instead of converting these massive routines wholesale, decompose them first:

1. Extract smaller, independent sections into separate kernels
2. Convert extracted kernels using index-based pattern
3. Leave complex cross-mesh sections (INTERPOLATED_BC, BACK_MESH) for later

**Example**:
- Extract SPECIFIED_TEMPERATURE case → separate kernel
- Extract CONVECTIVE_FLUX_BC case → separate kernel
- Leave INTERPOLATED_BC case in main routine (already uses OMESH, inherently cross-mesh)

This provides:
- ✅ Incremental progress
- ✅ Easier testing
- ✅ Reduced risk
- ✅ Better code organization

### Next Steps

**Option A**: Complete direct-access conversion of full routines (4-6 hours)
**Option B**: Decompose into smaller kernels first, then convert (2-3 hours per case)
**Option C**: Document current progress, focus on other parallelization targets

**Recommendation**: Option B for maintainability and reduced risk
