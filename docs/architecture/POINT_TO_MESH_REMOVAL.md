# POINT_TO_MESH Removal: Per-Module Difficulty Report

## The Problem

`POINT_TO_MESH(NM)` (mesh.f90:356-486) sets ~200 module-level pointer aliases
(`U => MESHES(NM)%U`, etc.). This is the single biggest barrier to intra-node
parallelism: if two threads call `POINT_TO_MESH` for different meshes, they
corrupt each other's pointers.

**Total occurrences: 171 calls across 21 files.**

## Summary Table

| Module | File | Lines | PTM Calls | Time-Step Critical? | Difficulty | Effort |
|--------|------|------:|----------:|:-------------------:|:----------:|:------:|
| MASS | mass.f90 | 57 | 0 | Yes | DONE | -- |
| DIVG | divg.f90 | 76 | 0 | Yes | DONE | -- |
| WALL_ROUTINES | wall.f90 | 3,239 | 2 | Yes | EASY | 1 day |
| SOOT | soot.f90 | 599 | 4 | Yes | EASY | 1 day |
| HVAC | hvac.f90 | 5,144 | 1 | Yes (rank 0) | EASY | 0.5 day |
| TURB | turb.f90 | 2,681 | 9 | Mostly init/test | EASY | 1 day |
| FIRE | fire.f90 | 1,153 | 9 | Yes | MODERATE | 2 days |
| VELO | velo.f90 | 2,109 | 9 | Yes | MODERATE | 3 days |
| PART | part.f90 | 4,850 | 4 | Yes | HARD | 4 days |
| VEGE | vege.f90 | 1,449 | 4 | Conditional | MODERATE | 2 days |
| PRES | pres.f90 | 5,205 | 30 | Yes | HARD | 5 days |
| RAD | radi.f90 | 4,838 | 1 | Yes | HARD | 4 days |
| DUMP | dump.f90 | 11,629 | 5 | Yes (output) | MODERATE | 2 days |
| INIT | init.f90 | 5,510 | 7 | Init only | LOW PRIO | -- |
| READ_INPUT | read.f90 | 17,013 | 11 | Init only | LOW PRIO | -- |
| COMPLEX_GEOMETRY | geom.f90 | 27,684 | 24 | Init (mostly) | VERY HARD | 10+ days |
| CC_INIT | ccib_init.f90 | 6,883 | 24 | Init only | LOW PRIO | -- |
| CC_VELOCITY | ccib_velocity.f90 | 4,648 | 9 | Yes (if CC_IBM) | HARD | 4 days |
| CC_DENSITY | ccib_density.f90 | 2,099 | 11 | Yes (if CC_IBM) | MODERATE | 2 days |
| CC_DIVERGENCE | ccib_divergence.f90 | 3,654 | 2 | Yes (if CC_IBM) | EASY | 1 day |
| CC_PRESSURE | ccib_pressure.f90 | 2,824 | 2 | Yes (if CC_IBM) | EASY | 1 day |
| CC_SCALARS | ccib.f90 | 653 | 2 | Mixed | EASY | 0.5 day |
| CC_EXCHANGE | ccib_exchange.f90 | 1,217 | 1 | Yes (if CC_IBM) | EASY | 0.5 day |
| CC_VERIFICATION | ccib_verification.f90 | 548 | 4 | Test only | LOW PRIO | -- |

---

## Detailed Per-Module Analysis

### DONE: mass.f90 (0 calls), divg.f90 (0 calls)

Already fully refactored. Pure delegation to `*_KERNELS` with `MESHES(NM)`.

---

### EASY: wall.f90 (2 calls, 3,239 lines)

**Routines with PTM:**
| Routine | Line | Lines | Description |
|---------|------|------:|-------------|
| WALL_BC | 53 | 226 | Main wall boundary condition routine |
| TGA_ANALYSIS | 3143 | 94 | TGA test mode (not in time loop) |

**WALL_BC** is the main routine. After `POINT_TO_MESH`, it:
1. Sets pointer aliases based on PREDICTOR flag (`UU=>US` or `UU=>U`, etc.)
2. Loops over wall cells, delegating most work to kernels that already take `MESHES(NM)`
3. The non-kernel work is pointer alias setup and OpenMP loop orchestration

**Approach**: Replace `POINT_TO_MESH` with `M => MESHES(NM)`, resolve the
PREDICTOR/CORRECTOR pointer aliases using `M%` access, pass `M` to any
remaining internal subroutines. The 6 module-level pointers (PBAR_P, RHOP,
UU, VV, WW, ZZP) become local pointers or direct `M%` references.

**Blockers**: Module-level pointer aliases (`PBAR_P`, `RHOP`, `UU`, `VV`, `WW`, `ZZP`)
are used by contained subroutines. These need to be converted to local pointers
or passed as arguments.

---

### EASY: soot.f90 (4 calls, 599 lines)

**Routines with PTM:**
| Routine | Line | Lines | Description |
|---------|------|------:|-------------|
| SETTLING_VELOCITY | 40 | 137 | Gravitational settling |
| CALC_AGGLOMERATION | 277 | 114 | Soot agglomeration |
| SOOT_SURFACE_OXIDATION | 408 | 128 | Surface oxidation |
| DROPLET_SCRUBBING | 551 | 45 | Droplet scrubbing |

All four routines are small, self-contained, and use mesh pointers for basic
array access (RHO, ZZ, MU, TMP, WALL, BOUNDARY_COORD, etc.). No cross-mesh
access. No complex control flow.

**Approach**: Add `TYPE(MESH_TYPE), INTENT(INOUT) :: M` as first argument,
replace all bare variable names with `M%` prefix. Straightforward sed-like
replacement.

**Blockers**: Module-level SAVE variables (BIN_S, BIN_M, etc.) and pointer
aliases (WC, CFA, B1, B2, BC). These would need to become local or be moved
into the mesh type.

---

### EASY: hvac.f90 (1 call, 5,144 lines)

**Routines with PTM:**
| Routine | Line | Lines | Description |
|---------|------|------:|-------------|
| HVAC_BC_IN | 2293 | 218 | Set HVAC boundary conditions from mesh data |

Only one routine uses PTM. HVAC_CALC itself runs on rank 0 only and doesn't
use mesh pointers. HVAC_BC_IN reads wall data from the mesh and writes to
HVAC node arrays.

**Approach**: Pass `MESHES(NM)` to HVAC_BC_IN, replace pointer variables with
`M%` access. The contained subroutine INITIALIZE_HVAC also accesses mesh data.

**Blockers**: None significant. HVAC solver runs on rank 0 only, so it's not
a parallelization target anyway. But HVAC_BC_IN is called per-mesh and should
be refactored.

---

### EASY: turb.f90 (9 calls, 2,681 lines)

**Routines with PTM:**
| Routine | Line | Lines | Description | Critical Path? |
|---------|------|------:|-------------|:-:|
| INIT_TURB_ARRAYS | 32 | 59 | Allocate turbulence work arrays | Init |
| NS_ANALYTICAL_SOLUTION | 105 | 48 | Set analytical velocity/pressure | Init/test |
| COMPRESSION_WAVE | 183 | 98 | Initialize compression wave test | Init/test |
| TWOD_VORTEX_CERFACS | 295 | 29 | 2D vortex initialization | Init/test |
| TWOD_VORTEX_UMD | 340 | 22 | 2D vortex initialization | Init/test |
| TWOD_SOBOROT_UMD | 383 | 50 | 2D rotation initialization | Init/test |
| SYNTHETIC_TURBULENCE | 764 | 186 | Synthetic eddy method | Time step |
| SAAD_MMS_1 | 2198 | 14 | MMS test case | Init/test |
| SHUNN_MMS_3 | 2236 | 18 | MMS test case | Init/test |

Only `SYNTHETIC_TURBULENCE` is on the time-stepping critical path; the rest
are initialization/verification routines. All use straightforward field access.

**Approach**: Pass `MESHES(NM)`, replace pointer names with `M%`. The
verification routines are particularly simple.

**Blockers**: None. COMPUTE_VISCOSITY (the main time-step routine) already
delegates to TURB_KERNELS without PTM.

---

### MODERATE: fire.f90 (9 calls, 1,153 lines)

**Routines with PTM:**
| Routine | Line | Lines | Description |
|---------|------|------:|-------------|
| COMBUSTION_LOAD_BALANCED | 60, 79 | 51 | Top-level combustion dispatcher |
| COMBUSTION_GENERAL_LOAD_BALANCED | 120, 249, 310 | 249 | Cell discovery + serial chemistry |
| DISTRIBUTE_CELLS_ACCROSS_MPI_PROCESSES | 574 | 90 | MPI load balancing |
| GATHER_CELLS_FROM_MPI_PROCESSES | 703 | 40 | Gather results back |
| COMBUSTION_BC | 822 | 21 | Set combustion boundary conditions |
| CONDENSATION_EVAPORATION | 870 | 278 | Condensation/evaporation model |

**Complexity**: The combustion system does **cross-mesh load balancing** --
it gathers chemically active cells from all meshes, redistributes them across
MPI processes, solves chemistry, then scatters results back. This is inherently
a global operation.

**Approach**:
- `COMBUSTION_BC` and `CONDENSATION_EVAPORATION` are straightforward per-mesh
  routines: pass `MESHES(NM)`, replace pointer vars.
- `COMBUSTION_LOAD_BALANCED` is a global orchestrator that loops over meshes
  internally. The PTM calls inside the mesh loop can be replaced with
  `M => MESHES(NM)` local pointer, but the cross-mesh load balancing design
  means this can't simply become per-mesh parallel.
- `DISTRIBUTE/GATHER_CELLS` access multiple meshes' data for load balancing.

**Blockers**: Cross-mesh load balancing architecture. SAVE variables
(T_CHEM_ODE, COMBUSTION_INIT, cell distribution arrays). The load-balanced
combustion is fundamentally a global operation and will remain a
synchronization point.

---

### MODERATE: velo.f90 (9 calls, 2,109 lines)

**Routines with PTM:**
| Routine | Line | Lines | Description | Kernel Exists? |
|---------|------|------:|-------------|:-:|
| VISCOSITY_BC | 65 | 44 | Viscosity boundary conditions | No |
| VELOCITY_FLUX | 134 | 44 | Delegates to kernel + adds CC terms | Yes (partial) |
| VELOCITY_FLUX_CYLINDRICAL | 200 | 117 | Full cylindrical velocity flux | No |
| NO_FLUX | 340 | 199 | Enforce div-free at solid boundaries | No |
| VELOCITY_PREDICTOR | 569 | 47 | Delegates to kernel + CC projection | Yes |
| VELOCITY_CORRECTOR | 643 | 52 | Delegates to kernel + CC/diagnostics | Yes |
| VELOCITY_BC | 742 | 770 | Velocity boundary conditions (largest) | No |
| MATCH_VELOCITY | 1546 | 211 | Match velocity at mesh boundaries | No |
| MATCH_VELOCITY_FLUX | 1787 | 132 | Match momentum flux at mesh boundaries | No |

**Complexity breakdown**:
- **Thin wrappers** (VELOCITY_FLUX, VELOCITY_PREDICTOR, VELOCITY_CORRECTOR):
  Delegate to kernels with some pre/post PTM setup. Easy to fix.
- **Boundary routines** (VISCOSITY_BC, VELOCITY_BC, NO_FLUX): Significant
  per-mesh computation using mesh pointer variables. VELOCITY_BC alone is 770
  lines of boundary condition logic. These need full kernel extraction.
- **Cross-mesh routines** (MATCH_VELOCITY, MATCH_VELOCITY_FLUX): Access
  `OMESH(NOM)` neighbor data. These are inherently inter-mesh and should
  remain synchronization points, but still need PTM removal to avoid corrupting
  other threads' pointers during the (sequential) exchange phase.

**Approach**: For the 3 thin wrappers, remove PTM and use `M => MESHES(NM)`.
For VELOCITY_BC, NO_FLUX, VISCOSITY_BC: extract into velo_kernels.f90.
For MATCH_VELOCITY/MATCH_VELOCITY_FLUX: use local `M => MESHES(NM)` pointer.

**Blockers**: VELOCITY_BC is large (770 lines) and uses many mesh pointer
variables. MATCH_VELOCITY accesses cross-mesh data via OMESH and EXTERNAL_WALL.

---

### HARD: part.f90 (4 calls, 4,850 lines)

**Routines with PTM:**
| Routine | Line | Lines | Description |
|---------|------|------:|-------------|
| INSERT_ALL_PARTICLES | 156 | 1,644 | Insert new particles (contains 7 sub-routines) |
| MOVE_PARTICLES | 1843 | 1,611 | Move particles through velocity field (contains 7 sub-routines) |
| PARTICLE_MASS_ENERGY_TRANSFER | 3612 | 940 | Heat/mass exchange with gas |
| PARTICLE_MOMENTUM_TRANSFER | 4570 | 36 | Apply particle drag to gas |

**Complexity**: Large routines with many deeply-nested contained subroutines.
Each top-level routine calls PTM once, then all contained subroutines use the
module-level pointers. The contained subroutines access `U, V, W, RHO, ZZ,
TMP, MU, WALL, BOUNDARY_COORD, LAGRANGIAN_PARTICLE, CELL, CELL_INDEX`,
and various other mesh fields.

Additionally, **cross-mesh access** exists in multiple routines:
- `INSERT_ALL_PARTICLES`: accesses `MESHES(NOM)` (lines ~223-232) and
  `OMESH(NOM)` (line ~295) for inter-mesh particle transfer
- `PARTICLE_MASS_ENERGY_TRANSFER`: accesses `EXTERNAL_WALL` (line ~4023) and
  `OMESH(NOM)` (line ~4747) for boundary cell data from neighbor meshes

**Approach**: Add `M` parameter to each top-level routine and pass it down
to contained subroutines via host association (Fortran contained subroutines
can access the parent's locals). Replace ~200+ mesh pointer variable
references with `M%` prefix. Cross-mesh access sections will need special
handling — either passed as additional arguments or kept as explicit
`MESHES(NOM)` references.

**Blockers**: Large volume of code (4,195 lines of PTM-dependent code).
Cross-mesh access via `OMESH`, `EXTERNAL_WALL`, and `MESHES(NOM)` in
particle insertion and mass/energy transfer. The main challenges are
the sheer number of pointer variable substitutions and resolving the
cross-mesh data access patterns.

---

### MODERATE: vege.f90 (4 calls, 1,449 lines)

**Routines with PTM:**
| Routine | Line | Lines | Description | Critical Path? |
|---------|------|------:|-------------|:-:|
| INITIALIZE_LEVEL_SET_FIRESPREAD_1 | 45 | 109 | Init (allocate arrays) | Init |
| INITIALIZE_LEVEL_SET_FIRESPREAD_2 | 171 | 154 | Init (set initial values) | Init |
| LEVEL_SET_FIRESPREAD | 348 | 298 | Main fire spread computation | Yes (if enabled) |
| UPDATE_FIRE_SPREAD_OUTPUTS | 1385 | 62 | Update output quantities | Yes (if enabled) |

Only active when `LEVEL_SET_MODE > 0`. Uses module-level SAVE pointers
(PHI_LS_P, M, WC, CFA, etc.) and mutable state.

**Approach**: Replace module-level pointers with local `M => MESHES(NM)`,
pass M to contained subroutines.

**Blockers**: Module-level SAVE variables and pointers. Only matters when
level-set fire spread is enabled.

---

### HARD: pres.f90 (30 calls, 5,205 lines)

**Routines with PTM (grouped by subsystem):**

| Subsystem | Routines | PTM Calls | Lines | Critical Path? |
|-----------|----------|----------:|------:|:-:|
| **Core solver** | COMPUTE_VELOCITY_ERROR | 1 | 250 | Yes |
| | PRESSURE_SOLVER_CHECK_RESIDUALS_U | 1 | ~50 | Yes |
| **ULMAT solver** | ULMAT_SOLVER_SETUP, ULMAT_SOLVER | 3 | 275 | Init + Yes |
| | FINISH_ULMAT_SOLVER | 1 | 28 | Finalization |
| **GLMAT solver** | GLMAT_SOLVER | 5 | 242 | Yes (if GLMAT) |
| | GLMAT_SOLVER_SETUP | 2 | 143 | Init |
| | CHECK_UNSUPPORTED_MESH | 1 | 34 | Init |
| **Cross-mesh copy** | COPY_H_OMESH_TO_MESH | 2 | 169 | Yes (if GLMAT) |
| | COPY_HS_IN_CCVAR | 3 | 56 | Yes (if CC_IBM) |
| | COPY_CCVAR_IN_HS | 1 | 22 | Yes (if CC_IBM) |
| **Matrix assembly** | GET_BCS_H_MATRIX | 1 | 74 | Init |
| | GET_H_MATRIX | 1 | 239 | Init |
| | GET_MATRIXGRAPH_H_WHLDOM | 3 | 170+ | Init |
| | GET_H_REGFACES | 1 | ~80 | Init |
| | GET_MATRIX_INDEXES_H | 2 | ~80 | Init |
| | SET_CCVAR_CGSC_H | 1 | ~40 | Init |
| **Diagnostic** | WRITE_EWC_TYPE_DIAGNOSTIC | 1 | 46 | Diagnostic |

**Complexity**: pres.f90 has 3 distinct pressure solver strategies (FFT,
ULMAT, GLMAT). FFT is already kernel-extracted. ULMAT and GLMAT have deep
PTM usage. GLMAT_SOLVER is inherently cross-mesh (solves a single global
matrix across all meshes). COMPUTE_VELOCITY_ERROR accesses OMESH neighbor data.

**Approach**:
- **FFT path** (PRESSURE_SOLVER_COMPUTE_RHS, PRESSURE_SOLVER_FFT,
  PRESSURE_SOLVER_CHECK_RESIDUALS): Already kernel-extracted. No PTM.
- **COMPUTE_VELOCITY_ERROR**: Large, uses cross-mesh data via OMESH.
  Could be extracted to a kernel but needs OMESH access pattern resolved.
- **ULMAT_SOLVER**: Per-mesh solve. Extract to kernel pattern.
- **GLMAT_SOLVER**: Cross-mesh by design. Cannot be made per-mesh parallel.
  Will remain a synchronization point.
- **Matrix assembly**: Init-only, can be deprioritized.

**Blockers**: GLMAT_SOLVER is inherently global. COMPUTE_VELOCITY_ERROR
uses cross-mesh neighbor data. SAVE variables (CYL_FCT, ILO/IHI_CELL/FACE)
throughout.

---

### HARD: radi.f90 (1 call, 4,838 lines)

**Routines with PTM:**
| Routine | Line | Lines | Description |
|---------|------|------:|-------------|
| COMPUTE_RADIATION | 3438 | 1,210 | Complete radiation transport solver |

**Complexity**: Only 1 PTM call, but it gates access to the entire
RADIATION_FVM subroutine (~1,130 lines) and ADD_VOLUMETRIC_HEAT_SOURCE (~48
lines) which are contained subroutines. RADIATION_FVM is the full finite-volume
radiation transport solver that loops over angles, bands, and spatial cells.
It uses extensive mesh pointer variables: `RHO, TMP, KAPPA_GAS, QR, QR_W, MU,
ZZ, WALL, BOUNDARY_COORD, BOUNDARY_PROP1, BOUNDARY_RADIA, CFACE`, plus the
large spectral/angular radiation arrays.

The module also has multiple SAVE arrays at module level for spectral data
(CPLXREF_WATER, CPLXREF_FUEL, absorption coefficients, etc.) and uses
the RADCAL_VAR module with shared state.

**Approach**: Add `TYPE(MESH_TYPE), INTENT(INOUT) :: M` parameter. The
contained subroutines (RADIATION_FVM, ADD_VOLUMETRIC_HEAT_SOURCE) access the
parent's locals via host association. Replace all mesh pointer variables with
`M%` prefix. This is a large but mechanical substitution.

**Blockers**: Volume of code (~1,200 lines of mesh pointer usage). Module-level
SAVE arrays for spectral data. RADCAL_VAR shared state.

---

### MODERATE: dump.f90 (5 calls, 11,629 lines)

**Routines with PTM:**
| Routine | Line | Lines | Description | Critical Path? |
|---------|------|------:|-------------|:-:|
| UPDATE_GLOBAL_OUTPUTS | 74 | 9 | Delegate to sub-routines | Yes (every step) |
| DUMP_MESH_OUTPUTS | 99 | 119 | Dispatch output writes | Yes (when output due) |
| DUMP_RESTART | 3542 | 155 | Write restart file | Periodic |
| READ_RESTART | 3728 | 218 | Read restart file | Init only |
| WRITE_CFACES | 5531 | 6 | Write cut-face geometry | Periodic |

**Complexity**: UPDATE_GLOBAL_OUTPUTS and DUMP_MESH_OUTPUTS are called every
time step. They dispatch to many sub-routines (UPDATE_HRR, UPDATE_MASS,
DUMP_SLCF, DUMP_BNDF, DUMP_PART, etc.) that access mesh data through the
aliased pointers. However, output routines don't need to be parallelized
(they're not compute-bound), so the fix is just to prevent pointer corruption
during the parallel computation phase.

**Approach**: Replace PTM with `M => MESHES(NM)` local pointer. The contained
sub-routines use host association to access M.

**Blockers**: Many contained subroutines (~30+) that all rely on mesh pointer
aliases. Large mechanical refactoring but not conceptually difficult.

---

### LOW PRIORITY: init.f90 (7 calls), read.f90 (11 calls)

Init-only modules. Not on the time-stepping critical path. Can be deferred
indefinitely since initialization is sequential.

---

### LOW PRIORITY: geom.f90 (24 calls, 27,684 lines)

Mostly initialization (SET_CUTCELLS_3D, GET_EXT_INB_CUTFACES_TO_CFACE, etc.).
A few routines (EXCHANGE_CC_NOADVANCE_INFO, BLOCK_CC_SOLID_EXTWALLCELLS)
may be called during time-stepping when obstructions change.

**Rating**: VERY HARD due to sheer size (27K lines), 32 SAVE variables, and
deep interdependencies. Deprioritize unless CC_IBM time-stepping requires it.

---

### CC Modules During Time-Stepping

| Module | PTM Calls | Difficulty | Notes |
|--------|----------:|:----------:|-------|
| ccib_divergence.f90 | 2 | EASY | CC_DIVERGENCE_PART_1 + CC_CHECK_DIVERGENCE |
| ccib_pressure.f90 | 2 | EASY | GET_H_CUTFACES + (comment says "assume PTM done") |
| ccib.f90 | 2 | EASY | CC_RHO0W_INTERP + CC_H_INTERP |
| ccib_exchange.f90 | 1 | EASY | Inside FILL_GCCUTCELL_SPECIES |
| ccib_velocity.f90 | 9 | HARD | CC_RESTORE_UVW, CC_MATCH_VELOCITY*, CC_NO_FLUX, CC_COMPUTE_VELOCITY_ERROR, etc. Cross-mesh access. |
| ccib_density.f90 | 11 | MODERATE | Many routines, but per-mesh operations |
| ccib_init.f90 | 24 | LOW PRIO | Init only |
| ccib_verification.f90 | 4 | LOW PRIO | Test only |

---

## Recommended Removal Order

Based on impact (time-step critical path), difficulty, and dependencies:

### Phase 1: Quick Wins (5 days)
Removes PTM from the core per-mesh computation loop for non-CC_IBM cases.

1. **wall.f90** (2 calls) - EASY, on critical path
2. **soot.f90** (4 calls) - EASY, on critical path
3. **hvac.f90** (1 call) - EASY, on critical path
4. **turb.f90** (9 calls) - EASY, mostly init/test but SYNTHETIC_TURBULENCE is time-step

### Phase 2: Core Solvers (8 days)
Removes PTM from the main numerical solver routines.

5. **velo.f90** (9 calls) - MODERATE, extract VELOCITY_BC/NO_FLUX to kernels
6. **fire.f90** (9 calls) - MODERATE, per-mesh parts easy, load balancing stays global
7. **dump.f90** (5 calls) - MODERATE, mechanical but large

### Phase 3: Expensive Modules (10 days)
The largest and most complex modules.

8. **part.f90** (4 calls) - HARD, cross-mesh access + large volume of substitution
9. **radi.f90** (1 call) - HARD, single huge routine with deep pointer usage
10. **pres.f90** (30 calls) - HARD, multiple solver strategies, some cross-mesh

### Phase 4: CC_IBM (if needed) (8 days)

11. **ccib_divergence.f90** (2 calls) - EASY
12. **ccib_pressure.f90** (2 calls) - EASY
13. **ccib.f90** (2 calls) - EASY
14. **ccib_exchange.f90** (1 call) - EASY
15. **ccib_density.f90** (11 calls) - MODERATE
16. **ccib_velocity.f90** (9 calls) - HARD (cross-mesh)

### Defer: Init/Test Only
- init.f90, read.f90, geom.f90, ccib_init.f90, ccib_verification.f90

---

## Effort Summary

| Phase | Modules | PTM Calls Removed | Est. Effort |
|-------|---------|------------------:|:-----------:|
| Phase 1 (Quick wins) | 4 | 16 | 5 days |
| Phase 2 (Core solvers) | 3 | 23 | 8 days |
| Phase 3 (Expensive) | 3 | 35 | 10 days |
| Phase 4 (CC_IBM) | 6 | 27 | 8 days |
| Deferred (init/test) | 5 | 70 | -- |
| **Total critical path** | **16** | **101** | **~31 days** |

After Phase 2, the main time-stepping loop would be PTM-free for standard
(non-CC_IBM) simulations, enabling per-mesh parallel execution with Hedgehog.
