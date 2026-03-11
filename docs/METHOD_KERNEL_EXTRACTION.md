# Kernel Extraction Methodology

## Goal

Extract computation routines from FDS Fortran modules into dedicated `*_kernels.f90` modules that take `TYPE(MESH_TYPE)` as an explicit first argument instead of relying on global `MESH_POINTERS` / `POINT_TO_MESH`. This makes routines thread-safe and parallelizable.

## When to Apply

Apply this methodology when a Fortran module contains computation routines that:
- Are called per-mesh in a loop (`DO NM = ...`)
- Use `POINT_TO_MESH(NM)` to set up module-level pointer aliases
- Operate only on the local mesh (no `OMESH`, `MESHES()` array, `EXTERNAL_WALL` access)
- Contain meaningful computation (not pure orchestration or I/O)

## Completed Modules

| Parent Module | Kernel Module | Routines Extracted |
|---------------|---------------|--------------------|
| `velo.f90` | `velo_kernels.f90` | VELOCITY_PREDICTOR_KERNEL, VELOCITY_CORRECTOR_KERNEL, VELOCITY_FLUX_KERNEL, COMPUTE_VISCOSITY_KERNEL, CHECK_STABILITY_KERNEL, BAROCLINIC_CORRECTION_KERNEL |
| `divg.f90` | `divg_kernels.f90` | DIVERGENCE_PART_1_KERNEL, DIVERGENCE_PART_2_KERNEL, CHECK_DIVERGENCE_KERNEL |
| `mass.f90` | `mass_kernels.f90` | MASS_FINITE_DIFFERENCES_NEW_KERNEL, DENSITY_KERNEL |
| `turb.f90` | `turb_kernels.f90` | WALE_VISCOSITY, WALL_MODEL, TAU_WALL_IJ, TEST_FILTER_LOCAL, EX2G3D_KERNEL, FILL_EDGES_KERNEL, ... |
| `fire.f90` | `fire_kernels.f90` | COMBUSTION_MODEL, CHECK_REACTION, GET_FLAME_TEMPERATURE, ... |
| `wall.f90` | `wall_kernels.f90` | CALCULATE_RHO_D_F, CALC_DEPOSITION, PYROLYSIS |
| `pres.f90` | `pres_kernels.f90` | PRESSURE_SOLVER_COMPUTE_RHS, PRESSURE_SOLVER_FFT, PRESSURE_SOLVER_CHECK_RESIDUALS |
| `ccib_divergence.f90` | `ccib_divergence_kernels.f90` | 8 routines |
| `ccib_velocity.f90` | `ccib_velocity_kernels.f90` | 6 routines |
| `ccib_pressure.f90` | `ccib_pressure_kernels.f90` | 7 routines |

## Step-by-Step Procedure

### Step 1: Identify Candidate Routines

In the parent module (e.g., `velo.f90`), find routines that:

```
Good candidates:
  ✓ No POINT_TO_MESH call inside the routine body
  ✓ No OMESH / MESHES array / EXTERNAL_WALL access
  ✓ Operates only on local mesh data
  ✓ Contains meaningful computation loops
  ✓ Called from multiple places (high reuse value)

Poor candidates:
  ✗ Calls POINT_TO_MESH(NM) internally
  ✗ Reads OMESH(NOM)%... for ghost cell data
  ✗ Accumulates into global OUTPUT_DATA arrays (Q_DOT, M_DOT)
  ✗ Pure orchestration (just calls other routines + timing)
  ✗ Contains Fortran I/O (WRITE, PRINT, OPEN)
```

**Method**: Search for `POINT_TO_MESH`, `OMESH`, `EXTERNAL_WALL`, and `MESHES(` within each subroutine. Routines without these are candidates.

### Step 2: Create the Kernel Module File

Create `Source/<module>_kernels.f90` (or `Source/ccib/<module>_kernels.f90` for CCIB):

```fortran
!> \brief Thread-safe computation kernels extracted from <MODULE>.
!> These routines take TYPE(MESH_TYPE) as an argument instead of relying on
!> MESH_POINTERS / POINT_TO_MESH, decoupling computation from global state.

MODULE <MODULE>_KERNELS

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES, ONLY: MESH_TYPE
! Add USE TYPES if referencing REACTION, SPECIES_MIXTURE, etc.
! Add other USE only as needed by extracted routines

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE
PUBLIC :: KERNEL_ROUTINE_1, KERNEL_ROUTINE_2, ...

CONTAINS

! Extracted routines go here

END MODULE <MODULE>_KERNELS
```

**Critical imports**:
- `PRECISION_PARAMETERS` — for `EB`, `FB` types
- `GLOBAL_CONSTANTS` — for physical/numerical constants (does NOT import TYPES)
- `MESH_VARIABLES, ONLY: MESH_TYPE` — for the mesh type definition
- `TYPES` — only if routine references derived types like `REACTION`, `SPECIES_MIXTURE`

### Step 3: Transform Each Routine

For each candidate routine, apply these transformations:

#### 3a. Add mesh argument

```fortran
! Before:
SUBROUTINE COMPUTE_SOMETHING(T, DT, NM)
INTEGER, INTENT(IN) :: NM
...
CALL POINT_TO_MESH(NM)

! After:
SUBROUTINE COMPUTE_SOMETHING_KERNEL(M, T, DT, NM)
TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
INTEGER, INTENT(IN) :: NM  ! Keep NM if used for indexing
```

- Add `M` as the **first** positional argument
- Use `INTENT(INOUT)` if the routine modifies mesh state (most cases)
- Add `TARGET` if the routine creates pointers to mesh components
- Keep `NM` if it's used for global array indexing (e.g., `T_USED(NM)`, `DT_NEW(NM)`)

#### 3b. Replace mesh variable references

| Global Pointer | Explicit Mesh Reference |
|----------------|-------------------------|
| `RHO(I,J,K)` | `M%RHO(I,J,K)` |
| `U(I,J,K)` | `M%U(I,J,K)` |
| `IBAR` | `M%IBAR` |
| `WALL(IW)` | `M%WALL(IW)` |
| `CUT_CELL(ICC)%IJK` | `M%CUT_CELL(ICC)%IJK` |
| `RHOP => RHO` | `RHOP => M%RHO` |
| `RHOP => RHOS` | `RHOP => M%RHOS` |

**Pointer aliases**: Resolve `PREDICTOR` flag-based aliasing:
```fortran
! Before (in orchestration, after POINT_TO_MESH):
IF (PREDICTOR) THEN
   RHOP => RHO    ! Module-level pointer set by POINT_TO_MESH
ELSE
   RHOP => RHOS
ENDIF

! After (in kernel):
IF (PREDICTOR) THEN
   RHOP => M%RHO
ELSE
   RHOP => M%RHOS
ENDIF
```

#### 3c. Remove POINT_TO_MESH and USE MESH_POINTERS

Delete `CALL POINT_TO_MESH(NM)` and `USE MESH_POINTERS` from the kernel routine.

#### 3d. Keep global constants unchanged

`PREDICTOR`, `CORRECTOR`, `CC_IBM`, `SOLID_PHASE_ONLY`, `N_TRACKED_SPECIES`, `RSC_T`, `RPR_T`, etc. are global constants — they do NOT need `M%` prefix.

### Step 4: Update Callers

Every call site must pass the mesh explicitly:

```fortran
! In parent module (e.g., velo.f90) — typically has POINT_TO_MESH(NM) earlier:
! Before:
CALL COMPUTE_SOMETHING(T, DT, NM)
! After:
CALL COMPUTE_SOMETHING_KERNEL(MESHES(NM), T, DT, NM)

! In another kernel module — M already available as argument:
! Before:
CALL COMPUTE_SOMETHING(T, DT, NM)
! After:
CALL COMPUTE_SOMETHING_KERNEL(M, T, DT, NM)
```

### Step 5: Remove from Parent Module

Delete the extracted routine from the parent module. Add a `USE` import:

```fortran
! In parent module:
USE <MODULE>_KERNELS, ONLY: COMPUTE_SOMETHING_KERNEL
```

### Step 6: Add to CMakeLists.txt

Add the kernel file **before** the parent module in `FDS_FORTRAN_SOURCES`:

```cmake
# In CMakeLists.txt (FDS_FORTRAN_SOURCES list):
${CMAKE_SOURCE_DIR}/Source/<module>_kernels.f90    # NEW — must come before parent
${CMAKE_SOURCE_DIR}/Source/<module>.f90
```

The kernel module must compile before the parent because the parent `USE`s it.

### Step 7: Handle Line Length

Fortran has a 132-character line limit (enforced by `gfortran -Werror=line-truncation`). Adding `MESHES(NM)` or `M` as first argument may push lines over. Use continuation:

```fortran
! Too long:
CALL SET_EXIMRHOZZLIM_3D(MESHES(NM), RHO_ZZ_LSRC_FACE_X, RHO_ZZ_LSRC_FACE_Y, RHO_ZZ_LSRC_FACE_Z)

! Fixed:
CALL SET_EXIMRHOZZLIM_3D(MESHES(NM), &
   RHO_ZZ_LSRC_FACE_X, RHO_ZZ_LSRC_FACE_Y, &
   RHO_ZZ_LSRC_FACE_Z)
```

### Step 8: Build and Verify

```bash
# Build both targets
cd build_hh
cmake --build . --target fds -j$(nproc)
cmake --build . --target fds_hh -j$(nproc)

# Run test case
cd ../test_cases/run_1mesh
mpiexec --oversubscribe -n 1 ../../build_hh/fds dancing_eddies_1mesh_short.fds

# Compare against baseline
diff dancing_eddies_1mesh_short_devc.csv \
     ../archive/saved_results/orig_1mesh/dancing_eddies_1mesh_short_devc.csv
```

Results must be byte-identical. If not, a variable was missed in the `M%` transformation.

## Common Pitfalls

| Issue | Symptom | Fix |
|-------|---------|-----|
| Missing `USE TYPES` | Compile error on `REACTION`, `SPECIES_MIXTURE` | Add `USE TYPES` to kernel module header |
| Missing `M%` prefix | Wrong results or crash | Search routine for bare variable names that are mesh members |
| Line too long | `gfortran -Werror=line-truncation` | Break with `&` continuation |
| Cyclic dependency | Compile error | Kernel module must NOT `USE` parent module |
| Missing `TARGET` | Runtime segfault when creating pointers | Add `TARGET` to `TYPE(MESH_TYPE)` declaration |
| Wrong CMake order | Link error | Kernel file must appear before parent in FDS_FORTRAN_SOURCES |
| `#ifdef WITH_SUNDIALS` | Build fails on conditional compilation | Preserve `#ifdef` blocks when extracting |

## Automation Checklist

For each routine to extract:

1. [ ] Verify no `POINT_TO_MESH`, `OMESH`, `EXTERNAL_WALL`, `MESHES(` in routine body
2. [ ] Create/update `*_kernels.f90` with module header and imports
3. [ ] Copy routine, add `M` as first arg with `TYPE(MESH_TYPE), INTENT(INOUT), TARGET`
4. [ ] Replace all mesh-pointer variables with `M%` prefix
5. [ ] Remove `CALL POINT_TO_MESH` and `USE MESH_POINTERS`
6. [ ] Add routine to `PUBLIC` list in kernel module
7. [ ] Update all callers to pass `MESHES(NM)` or `M`
8. [ ] Remove routine from parent module, add `USE *_KERNELS, ONLY:`
9. [ ] Add to CMakeLists.txt before parent module
10. [ ] Build both `fds` and `fds_hh` targets
11. [ ] Verify byte-identical CSV output against baselines

## Advanced Pattern: Index-Based Access (No Pointers)

### Problem

When converting routines that use pointer variables to access ALLOCATABLE array elements, gfortran requires the TARGET attribute:

```fortran
! This fails without TARGET attribute on M%WALL:
TYPE(MESH_TYPE), INTENT(INOUT) :: M
TYPE(WALL_TYPE), POINTER :: WC
WC => M%WALL(WALL_INDEX)  ! Error: target is neither TARGET nor POINTER
```

### Solution: Use Integer Indices Instead

Instead of creating pointer variables, store the indices and use direct array access:

```fortran
! OLD approach (requires TARGET):
TYPE(WALL_TYPE), POINTER :: WC
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
TYPE(BOUNDARY_PROP2_TYPE), POINTER :: B2

WC => M%WALL(WALL_INDEX)
B1 => M%BOUNDARY_PROP1(WC%B1_INDEX)
B2 => M%BOUNDARY_PROP2(WC%B2_INDEX)

B1%HEAT_TRANS_COEF = 2.0_EB * B1%K_G * B1%RDN

! NEW approach (index-based, no TARGET needed):
INTEGER :: B1_INDEX, B2_INDEX

B1_INDEX = M%WALL(WALL_INDEX)%B1_INDEX
B2_INDEX = M%WALL(WALL_INDEX)%B2_INDEX

M%BOUNDARY_PROP1(B1_INDEX)%HEAT_TRANS_COEF = &
   2.0_EB * M%BOUNDARY_PROP1(B1_INDEX)%K_G * M%BOUNDARY_PROP1(B1_INDEX)%RDN
```

### Automated Replacement with awk

For routines with many pointer references, use awk to automate the replacement:

```bash
cat > /tmp/replace_pointers.awk << 'AWKEOF'
BEGIN { in_func = 0 }

/^SUBROUTINE YOUR_ROUTINE_NAME/ { in_func = 1 }

in_func {
    gsub(/B1%/, "M%BOUNDARY_PROP1(B1_INDEX)%")
    gsub(/B2%/, "M%BOUNDARY_PROP2(B2_INDEX)%")
    gsub(/BC%/, "M%BOUNDARY_COORD(BC_INDEX)%")
}

/^END SUBROUTINE YOUR_ROUTINE_NAME/ { in_func = 0 }

{ print }
AWKEOF

awk -f /tmp/replace_pointers.awk your_file.f90 > /tmp/new_file.f90
mv /tmp/new_file.f90 your_file.f90
```

### Watch for Line Length

Index-based access is more verbose and may exceed Fortran's 132-character limit:

```fortran
! Too long (>132 chars):
IF (ALLOCATED(M%BOUNDARY_PROP1(B1_INDEX)%M_DOT_G_PP_ACTUAL)) THEN

! Fixed with continuation:
IF (ALLOCATED(M%BOUNDARY_PROP1(B1_INDEX)%M_DOT_G_PP_ACTUAL)) &
   THEN
```

Find long lines:
```bash
awk 'length($0)>132 {print NR":"length($0)}' your_file.f90
```

### Examples

**HEAT_TRANSFER_COEFFICIENT** (func.f90) - Converted using this pattern
**CALC_HVAC_BC** (wall.f90) - Simpler case, didn't need pointer replacement

See: `docs/WALL_BC_CONVERSIONS_SUMMARY.md` for detailed examples.

### Benefits

- ✅ No need to modify MESH_TYPE definition (add TARGET)
- ✅ Works with existing Fortran compilers
- ✅ Thread-safe (no module-level state)
- ✅ Explicit mesh parameter enables parallelization

### Trade-offs

- More verbose code (`M%BOUNDARY_PROP1(B1_INDEX)%X` vs `B1%X`)
- Watch for line length issues
- Need to track multiple index variables

Despite verbosity, this is the **recommended approach** for converting routines with complex pointer usage.
