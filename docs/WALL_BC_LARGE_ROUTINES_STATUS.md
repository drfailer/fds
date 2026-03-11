# SURFACE_HEAT_TRANSFER and CALCULATE_ZZ_F Thread-Safe Conversion

## Status: ✅ COMPLETED

### Solution: Pointer-Based Approach

**Key insight**: Use `TYPE(MESH_TYPE), POINTER :: M` instead of `INTENT(INOUT)`.

Since MESHES is declared with TARGET attribute:
```fortran
TYPE (MESH_TYPE), SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: MESHES
```

We can use pointer assignment to elements of M:
```fortran
SUBROUTINE SURFACE_HEAT_TRANSFER(NM, PREDICTOR_FLAG, T, SF, BC, B1, ...)
  TYPE(MESH_TYPE), POINTER :: M
  REAL(EB), POINTER, DIMENSION(:,:,:) :: UU, VV, WW, RHOP
  REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP
  REAL(EB), POINTER, DIMENSION(:,:) :: PBAR_P
  
  M => MESHES(NM)  ! M now points to element with TARGET
  
  ! Conditional pointer setup based on predictor/corrector phase
  IF (PREDICTOR_FLAG) THEN
     RHOP => M%RHOS
     ZZP => M%ZZS
     UU => M%US
     VV => M%VS
     WW => M%WS
     PBAR_P => M%PBAR_S
  ELSE
     RHOP => M%RHO
     ZZP => M%ZZ
     UU => M%U
     VV => M%V
     WW => M%W
     PBAR_P => M%PBAR
  ENDIF
  
  ! Now use UU, VV, WW, RHOP, ZZP, PBAR_P as before
  UN = UU(II,JJ,KK)
  ...
END SUBROUTINE
```

### Why This Works

1. MESHES is declared with TARGET attribute at module level
2. M => MESHES(NM) makes M point to an element that has TARGET
3. Components of M (like M%US, M%RHO, etc.) can be targeted by pointers
4. No modification to MESH_TYPE definition needed
5. Thread-safe: no module-level state, all data passed explicitly

### Conversions Completed

#### ✅ SURFACE_HEAT_TRANSFER (379 lines)
- Signature: `(NM, PREDICTOR_FLAG, T, SF, BC, B1, WALL_INDEX, CFACE_INDEX, PARTICLE_INDEX)`
- Call sites updated: 3 (WALL cells, CFACE cells, particles)
- Array prefixing: TMP, RSUM, MU, DX, DY, DZ, etc. → M%
- Fixed: 2 line truncation errors (> 132 chars)
- Global arrays corrected: SURFACE, LAGRANGIAN_PARTICLE_CLASS (not M% members)

#### ✅ CALCULATE_ZZ_F (413 lines)
- Signature: `(NM, PREDICTOR_FLAG, T, DT, WALL_INDEX, CFACE_INDEX, PARTICLE_INDEX)`
- Call sites updated: 3
- Same pattern as SURFACE_HEAT_TRANSFER

### Testing

✅ **Build**: Successful  
✅ **Tests**: dancing_eddies_1mesh_short  
✅ **Results**: Byte-identical to baseline (both devc and hrr CSV files)

### Impact

**Total lines converted to thread-safe**: 792 lines (379 + 413)  
**Call sites updated**: 6  
**Module-level pointer dependencies eliminated**: PBAR_P, RHOP, UU, VV, WW, ZZP

Combined with previous conversions:
- CALC_HVAC_BC: 52 lines, 2 call sites
- HEAT_TRANSFER_COEFFICIENT: ~175 lines, 15+ call sites

**Grand total**: ~1019 lines of thread-safe code, 23+ call sites updated

### Lessons Learned

1. **Don't overcomplicate**: Original HEAT_TRANSFER_COEFFICIENT used `M => MESHES(NM)` pattern successfully
2. **Check TARGET attribute**: When pointer assignment fails, verify the target has TARGET or is itself a POINTER
3. **POINTER vs INTENT(INOUT)**: For derived types with ALLOCATABLE components, POINTER parameter allows more flexible access
4. **Line length matters**: Fortran 132-char limit requires continuation for verbose M% array access
5. **Global vs mesh arrays**: SURFACE, SPECIES_MIXTURE, LAGRANGIAN_PARTICLE_CLASS are global, not mesh members

### Next Steps for WALL_BC Parallelization

With all major callees now thread-safe, the path forward is clear:

**Phase 1**: Three-phase decomposition of main WALL_BC
1. ASSIGN_GHOST_VALUE (sequential, OMESH reads)
2. WALL_BC_PROCESS_CELLS_KERNEL (parallel, 90% of cells)
3. Cross-mesh finalization (INTERPOLATED_BC, BACK_MESH, CONSUME_MASS)

**Estimated effort**: 2-3 hours for extraction and hedgehog integration

**Expected benefit**: ~90% of WALL_BC wall cells can be processed in parallel across meshes
