# WALL_BC Decomposition Analysis

## Current Structure

`WALL_BC(T,DT,NM)` is a 239-line orchestration routine in `wall.f90` (lines 29-273) that applies thermal, species, and density boundary conditions to wall cells, cut-cell faces, and particles.

### Execution Flow

```
CALL POINT_TO_MESH(NM)                    [Sets module-level pointer aliases]
  ↓
Pointer setup (UU,VV,WW,RHOP,ZZP,PBAR_P)  [Based on PREDICTOR flag]
  ↓
WALL_CELL_LOOP_0 (lines 110-121)
  └─ ASSIGN_GHOST_VALUE (external cells)  [READS OMESH - CROSS-MESH]
  └─ NEAR_SURFACE_GAS_VARIABLES_KERNEL    [KERNEL - CELL-LOCAL]
  └─ HEAT_TRANSFER_COEFFICIENT            [CELL-LOCAL]
  ↓
THIN_WALL_CELL_LOOP_0 (lines 123-139)
  └─ NEAR_SURFACE_GAS_VARIABLES_KERNEL    [KERNEL - CELL-LOCAL]
  └─ HEAT_TRANSFER_COEFFICIENT            [CELL-LOCAL]
  ↓
WALL_CELL_LOOP (lines 143-186)
  └─ SURFACE_HEAT_TRANSFER or SOLID_HEAT_TRANSFER
      ├─ SURFACE_HEAT_TRANSFER (non-INTERPOLATED_BC) [MOSTLY CELL-LOCAL]
      ├─ SURFACE_HEAT_TRANSFER (INTERPOLATED_BC)     [READS OMESH - CROSS-MESH]
      └─ SOLID_HEAT_TRANSFER                         [READS BACK_MESH - CROSS-MESH]
  └─ CALCULATE_RHO_D_F                               [KERNEL - CELL-LOCAL]
  └─ WALL_MODEL                                      [CELL-LOCAL]
  └─ CALC_DEPOSITION                                 [KERNEL - CELL-LOCAL]
  └─ CALC_HVAC_BC                                    [CELL-LOCAL]
  └─ CALCULATE_ZZ_F
      ├─ Non-CONSUME_MASS parts                      [MOSTLY CELL-LOCAL]
      └─ CONSUME_MASS section                        [WRITES OMESH - CROSS-MESH]
  └─ CALCULATE_RHO_F_KERNEL                          [KERNEL - CELL-LOCAL]
  ↓
SOLID_HEAT_TRANSFER (thin walls)                     [READS BACK_MESH - CROSS-MESH]
  ↓
CFACE_LOOP (lines 198-237)
  └─ (Similar structure to WALL_CELL_LOOP)
  ↓
PARTICLE_LOOP (lines 241-268)
  └─ (Similar structure to WALL_CELL_LOOP)
```

## Cross-Mesh Dependencies (BLOCKER for parallelization)

### 1. ASSIGN_GHOST_VALUE (lines 276-395)
- **Location**: WALL_CELL_LOOP_0, called for external wall cells only
- **OMESH access**: Lines 297-330
  - Reads `OMESH(EWC%NOM)%RHOS/RHO` and `OMESH(EWC%NOM)%ZZS/ZZ`
  - Interpolates density/species from neighboring mesh to ghost cells
  - Writes to local `RHOP(BC%II,BC%JJ,BC%KK)`, `ZZP(...)`, `TMP(...)`
- **Function**: Sets ghost cell values at mesh boundaries by interpolating from adjacent mesh
- **Parallelizability**: **MUST BE SEQUENTIAL** — reads from other meshes

### 2. SURFACE_HEAT_TRANSFER — INTERPOLATED_BC case (lines 616-775)
- **Location**: WALL_CELL_LOOP, conditional on boundary type
- **OMESH access**: Lines 619-740
  - Reads `OMESH(EWC%NOM)%RHOS/RHO`, `OMESH(EWC%NOM)%ZZS/ZZ`, `OMESH(EWC%NOM)%MU`
  - Computes species diffusion fluxes across mesh boundaries
  - Heavy interpolation and averaging across fine/coarse mesh transitions
- **Function**: Applies thermal BC for INTERPOLATED_BOUNDARY (mesh-to-mesh coupling)
- **Parallelizability**: **INTERPOLATED_BC must be sequential**; other BC types are cell-local

### 3. CALCULATE_ZZ_F — CONSUME_MASS section (lines 1115-1137 in routine)
- **Location**: WALL_CELL_LOOP, conditional on CORRECTOR and consumable obstructions
- **OMESH access**: Lines 1127-1133
  - Writes to `OMESH(EWC%NOM)%REAL_SEND_PKG8` (mass loss data sent to neighboring mesh)
  - Reads `MESHES(EWC%NOM)%OBSTRUCTION(...)%MASS`
- **Function**: Tracks mass consumed from obstructions that span mesh boundaries
- **Parallelizability**: **CONSUME_MASS must be sequential** — writes to OMESH; rest is cell-local

### 4. SOLID_HEAT_TRANSFER — BACK_MESH (lines 69-74, 99-100 in routine)
- **Location**: WALL_CELL_LOOP, conditional on thermally-thick surfaces
- **OMESH access**: Lines 69-74 (and throughout routine)
  - Reads from `MESHES(BACK_MESH)%WALL(...)`, `MESHES(BACK_MESH)%BOUNDARY_PROP1(...)`
  - For back-to-back wall cells (thin partitions), reads boundary data from opposite mesh
- **Function**: Couples heat conduction through thin walls that span mesh boundaries
- **Parallelizability**: **BACK_MESH access must be sequential** — reads from other mesh's wall data

## Already-Extracted Kernels (CELL-LOCAL)

These routines have already been extracted to `wall_kernels.f90` and take `TYPE(MESH_TYPE)` as an explicit argument:

1. **NEAR_SURFACE_GAS_VARIABLES_KERNEL** (142 lines)
   - Computes gas-phase properties near wall surface
   - Helper routines: SCALAR_TO_POINT_K, GET_TRILINEAR_WEIGHTS_K
   - Used in: WALL_CELL_LOOP_0, THIN_WALL_CELL_LOOP_0, PARTICLE_LOOP

2. **CALCULATE_RHO_D_F** (in wall_kernels.f90)
   - Computes species diffusion coefficients at wall
   - Used in: WALL_CELL_LOOP, CFACE_LOOP

3. **CALCULATE_RHO_F_KERNEL** (60 lines)
   - Computes boundary density from species mass fractions
   - Used in: WALL_CELL_LOOP, CFACE_LOOP

4. **CALC_DEPOSITION** (in wall_kernels.f90)
   - Applies deposition BC for particles settling on walls
   - Used in: WALL_CELL_LOOP, CFACE_LOOP

5. **PYROLYSIS** (in wall_kernels.f90, called from SOLID_HEAT_TRANSFER)
   - 1-D pyrolysis model for solid decomposition
   - Cell-local when not using BACK_MESH

## Proposed Decomposition

### Design Principle

Separate WALL_BC into **three phases**:

1. **Sequential Pre-Processing** (OMESH reads)
2. **Parallel Kernel Execution** (cell-local computation)
3. **Sequential Post-Processing** (OMESH writes)

### Phase 1: WALL_BC_GHOST_INTERPOLATION (Sequential)

**Purpose**: Set up ghost cell values by reading from neighboring meshes.

**Processes**:
- Loop through all **EXTERNAL** wall cells across all meshes
- Call `ASSIGN_GHOST_VALUE` to interpolate density/species from OMESH
- Writes to ghost cells: `RHO(BC%II,BC%JJ,BC%KK)`, `ZZ(...)`, `TMP(...)`

**Why Sequential**: Reads `OMESH(NOM)%RHOS/RHO/ZZS/ZZ` from neighboring meshes.

**Execution**:
```fortran
DO NM = 1, NMESHES
  CALL POINT_TO_MESH(NM)
  DO IW = 1, N_EXTERNAL_WALL_CELLS
    IF (external wall) CALL ASSIGN_GHOST_VALUE(IW,BC,B1)
  ENDDO
ENDDO
```

**Parallelization**: Each mesh processed sequentially (cannot parallelize across meshes due to OMESH reads).

---

### Phase 2: WALL_BC_PROCESS_CELLS_KERNEL (Parallel)

**Purpose**: Apply cell-local boundary conditions (thermal, species, density) to all wall cells.

**Signature**:
```fortran
SUBROUTINE WALL_BC_PROCESS_CELLS_KERNEL(M, T, DT, &
  CALL_HT_1D, PREDICTOR_FLAG, WALL_COUNTER, BC_CLOCK)
  TYPE(MESH_TYPE), INTENT(INOUT) :: M
  REAL(EB), INTENT(IN) :: T, DT
  LOGICAL, INTENT(IN) :: CALL_HT_1D, PREDICTOR_FLAG
  INTEGER, INTENT(IN) :: WALL_COUNTER
  REAL(EB), INTENT(IN) :: BC_CLOCK
```

**Processes**:

1. **Near-Surface Gas Variables** (all cells):
   ```fortran
   DO IW = 1, N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS
     ! Skip ASSIGN_GHOST_VALUE (done in Phase 1)
     CALL NEAR_SURFACE_GAS_VARIABLES_KERNEL(M, T, SF, BC, B1, WALL_INDEX=IW)
     IF (CALL_HT_1D) B1%HEAT_TRANS_COEF = HEAT_TRANSFER_COEFFICIENT(...)
   ENDDO
   ```

2. **Thermal Boundary Conditions** (non-INTERPOLATED, non-BACK_MESH):
   ```fortran
   DO IW = 1, N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS
     IF (WC%BOUNDARY_TYPE /= INTERPOLATED_BOUNDARY) THEN
       IF (.NOT.SF%THERMAL_BC_INDEX==THERMALLY_THICK) THEN
         CALL SURFACE_HEAT_TRANSFER_LOCAL(...)  ! Excludes INTERPOLATED_BC path
       ELSEIF (CALL_HT_1D .AND. .NOT. has_BACK_MESH) THEN
         CALL SOLID_HEAT_TRANSFER_LOCAL(...)    ! Excludes BACK_MESH path
       ENDIF
     ENDIF
   ENDDO
   ```

3. **Species and Density BC** (non-CONSUME_MASS):
   ```fortran
   DO IW = 1, N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS
     CALL CALCULATE_RHO_D_F(M, B1, BC, WALL_INDEX=IW)
     CALL WALL_MODEL(...)
     CALL CALC_DEPOSITION(M, DT, BC, B1, B2, WALL_INDEX=IW)
     CALL CALC_HVAC_BC(BC, B1, SF)
     CALL CALCULATE_ZZ_F_LOCAL(...)  ! Excludes CONSUME_MASS section
     CALL CALCULATE_RHO_F_KERNEL(M, BC, B1, WALL_INDEX=IW)
   ENDDO
   ```

4. **CFACEs and Particles** (similar structure):
   ```fortran
   DO ICF = ..., N_INTERNAL_CFACE_CELLS
     ! (Same as wall cells, minus ASSIGN_GHOST_VALUE)
   ENDDO
   DO IP = 1, NLP
     ! (Same as wall cells, minus cross-mesh operations)
   ENDDO
   ```

**What's Excluded** (deferred to Phases 1 & 3):
- ASSIGN_GHOST_VALUE → Phase 1
- SURFACE_HEAT_TRANSFER (INTERPOLATED_BC case) → Phase 3
- SOLID_HEAT_TRANSFER (BACK_MESH case) → Phase 3
- CALCULATE_ZZ_F (CONSUME_MASS section) → Phase 3

**Why Parallel**: All operations are cell-local — no OMESH reads/writes.

**Thread Safety**:
- Requires `TYPE(MESH_TYPE)` argument (no `POINT_TO_MESH`)
- Each thread processes a different mesh
- No shared state between threads

---

### Phase 3: WALL_BC_CROSS_MESH_FINALIZE (Sequential)

**Purpose**: Process wall cells with cross-mesh dependencies (INTERPOLATED_BC, BACK_MESH, CONSUME_MASS).

**Processes**:

1. **INTERPOLATED_BC Thermal Coupling**:
   ```fortran
   DO NM = 1, NMESHES
     CALL POINT_TO_MESH(NM)
     DO IW = 1, N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS
       IF (WC%BOUNDARY_TYPE == INTERPOLATED_BOUNDARY) THEN
         CALL SURFACE_HEAT_TRANSFER(NM, T, SF, BC, B1, WALL_INDEX=IW)
       ENDIF
     ENDDO
   ENDDO
   ```

2. **BACK_MESH Solid Heat Transfer**:
   ```fortran
   DO NM = 1, NMESHES
     CALL POINT_TO_MESH(NM)
     DO IW = 1, N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS
       IF (CALL_HT_1D .AND. has_BACK_MESH) THEN
         CALL SOLID_HEAT_TRANSFER(NM, T, DT_BC, WALL_INDEX=IW)
       ENDIF
     ENDDO
     DO ITW = 1, N_THIN_WALL_CELLS
       CALL SOLID_HEAT_TRANSFER(NM, T, 3*DT_BC, THIN_WALL_INDEX=ITW)
     ENDDO
   ENDDO
   ```

3. **CONSUME_MASS** (writes to OMESH):
   ```fortran
   DO NM = 1, NMESHES
     CALL POINT_TO_MESH(NM)
     DO IW = 1, N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS
       IF (CORRECTOR .AND. consumable_obstruction) THEN
         CALL CALCULATE_ZZ_F_CONSUME_MASS(NM, T, DT, WALL_INDEX=IW)
       ENDIF
     ENDDO
   ENDDO
   ```

**Why Sequential**: Reads/writes OMESH and MESHES(other_NM) data.

**Execution**: Each mesh processed sequentially.

---

## Parallelization Potential

### Phase 2 Kernel Coverage

**Total wall cells** (typical FDS case):
- EXTERNAL walls: ~40% of total walls
  - INTERPOLATED_BC: ~10-15% of external walls (mesh boundaries with mismatched grids)
  - Standard BC: ~85-90% of external walls
- INTERNAL walls: ~60% of total walls (all standard BC)

**Parallel coverage**: ~90-95% of wall cells can be processed in Phase 2 kernel.

**BACK_MESH cells**: Rare (thin partitions spanning meshes). Typically < 1% of walls.

**CONSUME_MASS cells**: Only external walls with consumable obstructions. Typically < 5% of walls.

### Speedup Estimate

Assuming:
- 4 meshes, 4 kernel threads
- 90% of wall processing in Phase 2 (parallel)
- 10% in Phases 1+3 (sequential)

**Amdahl's Law**:
- Sequential fraction: 10%
- Max speedup: 1 / (0.10 + 0.90/4) ≈ **3.1×**

Combined with existing parallelized sub-graphs (velocity, divergence, density), WALL_BC decomposition could reduce the current 49% sequential bottleneck to **~30%**, raising the overall max speedup from **2.05×** to **~2.9×**.

---

## Summary

| Phase | Operations | OMESH Access | Parallelizable |
|-------|-----------|--------------|----------------|
| **1. Ghost Interpolation** | ASSIGN_GHOST_VALUE | ✓ Reads OMESH | ❌ Sequential |
| **2. Process Cells** | Near-surface gas vars, thermal BC (non-INTERPOLATED), species BC (non-CONSUME_MASS) | ✗ Cell-local | ✅ Parallel kernel |
| **3. Cross-Mesh Finalize** | INTERPOLATED_BC, BACK_MESH, CONSUME_MASS | ✓ Reads/writes OMESH | ❌ Sequential |

**Key Insight**: 90-95% of wall cell processing is cell-local and can be parallelized. The 5-10% requiring cross-mesh access must remain sequential, but this is acceptable given the large coverage of the parallel kernel.
