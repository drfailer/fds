# WALL_BC "Quick Wins" - Thread Safety Status

## Question: Are the quick wins already kernels?

**Answer**: Mixed — one is almost there, one needs conversion, one is already thread-safe.

## Detailed Analysis

### 1. HEAT_TRANSFER_COEFFICIENT (func.f90) — ⚠️ 95% KERNEL

**Current signature:**
```fortran
REAL(EB) FUNCTION HEAT_TRANSFER_COEFFICIENT(NM, T, DELTA_N_TMP, SF, ...)
  INTEGER, INTENT(IN) :: NM
  ...
  TYPE(MESH_TYPE), POINTER :: M

  M => MESHES(NM)  ! ← Still relies on MESHES array
```

**Thread safety status**: NOT SAFE
- Takes mesh number `NM` as argument
- Uses `M => MESHES(NM)` internally
- All other access is via `M%...` (good!)

**What needs to change**:
```fortran
REAL(EB) FUNCTION HEAT_TRANSFER_COEFFICIENT(M, T, DELTA_N_TMP, SF, ...)
  TYPE(MESH_TYPE), INTENT(IN) :: M  ! ← Pass mesh explicitly
  ! Remove: M => MESHES(NM)
```

**Effort**: 20 minutes (search/replace `M =>` assignment, update call sites)

---

### 2. CALC_HVAC_BC (wall.f90) — ❌ NOT A KERNEL

**Current signature:**
```fortran
SUBROUTINE CALC_HVAC_BC(BC, B1, SF)
```

**Thread safety status**: NOT SAFE
- Uses module-level pointer `PBAR_P` (set by WALL_BC based on PREDICTOR flag)
- Line: `B1%RHO_F = PBAR_P(BC%KK, B1%PRESSURE_ZONE) / (RSUM_F * B1%TMP_G)`

**What needs to change**:
```fortran
SUBROUTINE CALC_HVAC_BC(M, PREDICTOR_FLAG, BC, B1, SF)
  TYPE(MESH_TYPE), INTENT(IN) :: M
  LOGICAL, INTENT(IN) :: PREDICTOR_FLAG
  ...
  ! Replace PBAR_P with:
  IF (PREDICTOR_FLAG) THEN
    B1%RHO_F = M%PBAR_S(BC%KK, B1%PRESSURE_ZONE) / (RSUM_F * B1%TMP_G)
  ELSE
    B1%RHO_F = M%PBAR(BC%KK, B1%PRESSURE_ZONE) / (RSUM_F * B1%TMP_G)
  ENDIF
```

**Effort**: 15 minutes (add arguments, conditional replacement, update call sites)

---

### 3. DEPOSIT_PARTICLE_MASS (wall.f90) — ✅ ALREADY A KERNEL

**Current signature:**
```fortran
SUBROUTINE DEPOSIT_PARTICLE_MASS(LP, LPC)
```

**Thread safety status**: SAFE
- Does NOT use module-level pointers (PBAR_P, RHOP, UU, VV, WW, ZZP)
- Accesses global arrays (SURFACE, BOUNDARY_COORD, BOUNDARY_PROP1) via indices from LP/LPC
- These global arrays are read-only during wall BC processing

**What needs to change**: NOTHING for thread safety

**Optional improvement**: Move to `wall_kernels.f90` for better organization:
```fortran
SUBROUTINE DEPOSIT_PARTICLE_MASS_KERNEL(M, LP, LPC)
  TYPE(MESH_TYPE), INTENT(IN) :: M
  ! Replace SURFACE with M%SURFACE (if needed)
  ! Replace BOUNDARY_COORD with M%BOUNDARY_COORD
  ! Replace BOUNDARY_PROP1 with M%BOUNDARY_PROP1
```

**Effort**: 15 minutes (move to wall_kernels.f90, add M argument for consistency)

---

## Summary

| Routine | Thread Safe? | In Kernel Module? | Effort to Convert |
|---------|--------------|-------------------|-------------------|
| HEAT_TRANSFER_COEFFICIENT | ⚠️ Almost (95%) | No (in func.f90) | 20 min |
| CALC_HVAC_BC | ❌ No (uses PBAR_P) | No (in wall.f90) | 15 min |
| DEPOSIT_PARTICLE_MASS | ✅ Yes | No (in wall.f90) | 0 min (optional: move to wall_kernels) |

**Total effort for thread safety**: 35 minutes

**Total effort including organization**: 50 minutes (if we also move DEPOSIT_PARTICLE_MASS to wall_kernels.f90)

## Note on Global Arrays

Routines that access global arrays like `SURFACE`, `BOUNDARY_COORD`, `BOUNDARY_PROP1` via mesh-independent indices are still thread-safe AS LONG AS:
1. The arrays are read-only during the parallel section
2. Different threads work on different indices (no write conflicts)

`DEPOSIT_PARTICLE_MASS` satisfies both conditions — it reads from these arrays using particle-specific indices (LP%BC_INDEX, LP%B1_INDEX) and writes only to particle-local data.

## Recommendation

Convert all three for consistency and to enable future parallelization:
1. **HEAT_TRANSFER_COEFFICIENT** - Critical dependency for SURFACE_HEAT_TRANSFER
2. **CALC_HVAC_BC** - Simple replacement
3. **DEPOSIT_PARTICLE_MASS** - Move to wall_kernels.f90 for organization

This sets up a clean foundation for the harder conversions (SURFACE_HEAT_TRANSFER, CALCULATE_ZZ_F).
