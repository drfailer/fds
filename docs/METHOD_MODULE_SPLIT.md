# Module Split Methodology

## Goal

Decompose large monolithic Fortran modules into smaller, focused modules organized by functional area. This reduces compilation coupling, improves maintainability, and enables targeted kernel extraction.

## When to Apply

Apply when a module:
- Exceeds ~5,000 lines and contains routines spanning multiple functional areas
- Has routines with different parallelization characteristics (some thread-safe, some not)
- Would benefit from independent compilation of subsets

## Completed Splits

| Original Module | Lines | Result | Modules Created |
|-----------------|-------|--------|-----------------|
| `ccib.f90` (CC_SCALARS) | 23,000 | 8 modules | ccib_data, ccib_init, ccib_exchange, ccib_divergence, ccib_density, ccib_velocity, ccib_pressure, ccib_verification |

## Step-by-Step Procedure

### Step 1: Survey the Module

Catalog all routines and classify them by functional area:

```
Module: CC_SCALARS (23K lines)
├── Data declarations (SAVE variables, parameters, pointers)
├── Init/Setup (18 routines, ~6800 lines)
├── MPI Exchange (1 routine, ~1200 lines)
├── Divergence (11 routines, ~4060 lines)
├── Density/Species (10 routines, ~2062 lines)
├── Velocity (24 routines, ~5109 lines)
├── Pressure (solver routines, ~2955 lines)
├── Verification (test routines, ~548 lines)
└── Shared helpers (~664 lines)
```

**Method**: grep for `SUBROUTINE` and `FUNCTION` declarations, then group by prefix/purpose.

### Step 2: Plan the Split Order

Split in dependency order — modules with fewer dependencies first:

1. **Data module** first (shared declarations, no routines)
2. **Leaf modules** next (routines that don't call other routines in the module)
3. **Dependent modules** last (routines that call routines in other new modules)

The goal is to maintain an **acyclic** dependency graph. If two modules would create a circular dependency, keep the shared routines in a common helper module.

### Step 3: Extract the Data Module

Create a dedicated data module for module-level declarations:

```fortran
MODULE CC_SCALARS_DATA

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
! ... other USE statements

IMPLICIT NONE (TYPE,EXTERNAL)

! Move all module-level declarations here:
! - SAVE variables
! - Parameters
! - Convenience pointers
! - Derived type definitions local to this module

INTEGER, SAVE :: MODULE_VARIABLE_1
REAL(EB), SAVE :: MODULE_VARIABLE_2
! ...

END MODULE CC_SCALARS_DATA
```

**All new modules** will `USE CC_SCALARS_DATA` instead of relying on host-association from the parent module.

### Step 4: Extract Each Functional Module

For each functional area, create a new module file:

#### 4a. Create the new file

```fortran
MODULE CC_DIVERGENCE

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES
USE TYPES
USE CC_SCALARS_DATA          ! Shared data from parent
! Add other USE as needed by extracted routines

IMPLICIT NONE (TYPE,EXTERNAL)

PRIVATE
PUBLIC :: ROUTINE_1, ROUTINE_2, ...

CONTAINS

! Paste extracted routines here

END MODULE CC_DIVERGENCE
```

#### 4b. Move routines from parent to new module

- Cut each routine from the parent module
- Paste into the new module's `CONTAINS` section
- Adjust `USE` statements in the new module header to cover all dependencies

#### 4c. Update the parent module

In the original module, remove the extracted routines and add:
```fortran
USE CC_DIVERGENCE, ONLY: ROUTINE_1, ROUTINE_2, ...
```

If other modules directly call the extracted routines, update their `USE` statements too.

#### 4d. Handle cross-module calls

If an extracted routine calls another routine still in the parent:
- Option A: Extract both together into the same new module
- Option B: Make the called routine `PUBLIC` in the parent, and add `USE parent_module, ONLY: called_routine` in the new module
- Option C: Move the called routine to a shared helper module

**Avoid circular USE**: If module A needs routines from module B and vice versa, move the shared routines into a third module (the data or helper module).

### Step 5: Update Callers Throughout the Codebase

Search for all files that `USE` the original module or call the extracted routines:

```bash
# Find all callers
grep -rn "USE CC_SCALARS" Source/ --include="*.f90"
grep -rn "CALL ROUTINE_1" Source/ --include="*.f90"
```

Update each caller to import from the new module instead:

```fortran
! Before:
USE CC_SCALARS, ONLY: CC_DIVERGENCE_PART_1

! After:
USE CC_DIVERGENCE, ONLY: CC_DIVERGENCE_PART_1
```

### Step 6: Add to CMakeLists.txt

Add new modules in dependency order:

```cmake
# In FDS_FORTRAN_SOURCES:
${CMAKE_SOURCE_DIR}/Source/ccib/ccib_data.f90         # Data module first
${CMAKE_SOURCE_DIR}/Source/ccib/ccib.f90               # Shared helpers
${CMAKE_SOURCE_DIR}/Source/ccib/ccib_init.f90          # Leaf modules
${CMAKE_SOURCE_DIR}/Source/ccib/ccib_exchange.f90
${CMAKE_SOURCE_DIR}/Source/ccib/ccib_divergence.f90
${CMAKE_SOURCE_DIR}/Source/ccib/ccib_density.f90
${CMAKE_SOURCE_DIR}/Source/ccib/ccib_velocity.f90
${CMAKE_SOURCE_DIR}/Source/ccib/ccib_pressure.f90
${CMAKE_SOURCE_DIR}/Source/ccib/ccib_verification.f90
```

### Step 7: Build and Verify

```bash
cmake --build . --target fds -j$(nproc)
cmake --build . --target fds_hh -j$(nproc)

# Run test and compare against baseline
cd ../test_cases/run_1mesh
mpiexec --oversubscribe -n 1 ../../build_hh/fds dancing_eddies_1mesh_short.fds
diff dancing_eddies_1mesh_short_devc.csv ../archive/saved_results/orig_1mesh/dancing_eddies_1mesh_short_devc.csv
```

Results must be byte-identical — a module split should not change any behavior.

## Incremental Strategy

Split one functional area per commit. This makes it easy to bisect if something breaks:

```
Commit 1: Move file + extract data module
Commit 2: Extract init/setup routines
Commit 3: Extract MPI exchange routines
Commit 4: Extract divergence routines
Commit 5: Extract density/species routines
Commit 6: Extract velocity routines
Commit 7: Extract pressure + verification routines
```

Each commit should build and pass tests independently.

## Resulting Structure

After splitting `ccib.f90`:

```
Source/ccib/
├── ccib_data.f90              (CC_SCALARS_DATA — shared declarations)
├── ccib.f90                   (CC_SCALARS — shared helpers, ~664 lines)
├── ccib_init.f90              (CC_INIT — initialization)
├── ccib_exchange.f90          (CC_EXCHANGE — MPI exchange)
├── ccib_divergence.f90        (CC_DIVERGENCE — divergence)
├── ccib_density.f90           (CC_DENSITY_MOD — density/species)
├── ccib_velocity.f90          (CC_VELOCITY — velocity)
├── ccib_pressure.f90          (CC_PRESSURE — pressure solver)
└── ccib_verification.f90      (CC_VERIFICATION — test cases)
```

Dependency graph (acyclic):
```
CC_SCALARS_DATA ←── all modules
CC_SCALARS (helpers) ←── divergence, velocity, pressure
CC_INIT ──→ CC_VELOCITY (one-way)
```

## Relationship to Kernel Extraction

Module split is a **prerequisite** for kernel extraction in large modules:

1. **Split** the monolithic module into functional areas (this methodology)
2. **Extract kernels** from each functional module (see METHOD_KERNEL_EXTRACTION.md)

For smaller modules (velo.f90, divg.f90, mass.f90), kernel extraction can be done directly without a prior module split.

## Common Pitfalls

| Issue | Symptom | Fix |
|-------|---------|-----|
| Circular USE | Compile error: "module not found" | Move shared routines to helper/data module |
| Missing PUBLIC | Compile error: "not accessible" | Add routine to `PUBLIC` list in new module |
| Wrong file order in CMake | Compile error: "module file not found" | Data module first, then helpers, then dependents |
| Forgotten caller update | Link error or wrong module reference | Search codebase for all callers and update USE statements |
| Private helper routines | Compile error when called from new module | Move helper alongside its caller, or make PUBLIC |

## Automation Checklist

For each functional area to extract:

1. [ ] Identify all routines belonging to the functional area
2. [ ] Check for cross-references with other functional areas
3. [ ] Create new module file with proper header and USE statements
4. [ ] Move routines from parent to new module
5. [ ] Update parent: remove routines, add USE import
6. [ ] Update all callers across the codebase
7. [ ] Add to CMakeLists.txt in correct dependency order
8. [ ] Build both `fds` and `fds_hh` targets
9. [ ] Verify byte-identical CSV output
10. [ ] Commit with descriptive message
