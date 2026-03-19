# FDS Pipelining Parallelism Analysis

## Motivation

The current Hedgehog graph parallelizes FDS over meshes (each mesh processed on a separate thread) and within meshes via K-block decomposition. Both approaches have hit diminishing returns:

- **Mesh-level parallelism** is limited by the number of meshes and sequential barriers (mesh exchanges, global reductions).
- **K-block decomposition** is limited to routines with clean K-separability; DIV_PART_1 was found not viable due to interleaved wall loops and face-value data races.

**Pipelining** is a third axis of parallelism: running independent computation stages *concurrently within each timestep*. Instead of executing VELOCITY_FLUX, then WALL_BC, then DIVERGENCE sequentially, we identify stages with non-overlapping data dependencies and run them in parallel.

## Diagrams

| File | Description |
|------|-------------|
| [fds_dataflow.dot](fds_dataflow.dot) / [.svg](fds_dataflow.svg) | Complete data-flow dependency graph for all predictor and corrector routines |
| [fds_pipeline_opportunities.dot](fds_pipeline_opportunities.dot) / [.svg](fds_pipeline_opportunities.svg) | Focused view of the identified pipelining opportunities with dependency proof |

Regenerate SVGs with: `dot -Tsvg fds_dataflow.dot -o fds_dataflow.svg`

## Data Categories

The analysis groups MESH_TYPE arrays into 15 logical categories:

| Category | Arrays | Description |
|----------|--------|-------------|
| vel | U, V, W | Current velocity |
| vel_s | US, VS, WS | Estimated/predicted velocity |
| flux | FVX, FVY, FVZ (+ baroclinic) | Momentum flux terms |
| drag | FVX_D, FVY_D, FVZ_D | Particle drag forces on gas |
| pres | H | Current pressure head |
| pres_s | HS | Estimated pressure head |
| div | D | Current divergence |
| div_s | DS | Estimated divergence |
| dddt | DDDT | dD/dt time derivative |
| species / species_s | ZZ / ZZS | Species mass fractions |
| rho / rho_s | RHO / RHOS | Density |
| temp | TMP | Temperature |
| visc | MU, KRES | Viscosity, kinetic energy |
| energy | Q, M_DOT_PPP | Energy source, mass production |
| qr | QR | Radiation source term |
| pbar | PBAR_S, D_PBAR_DT_S | Pressure zone data |
| wall | WALL(:), BCs | Wall/boundary cell data |
| part | LAGRANGIAN_PARTICLE | Lagrangian particles |
| diff | DEL_RHO_D_DEL_Z | Species diffusion fluxes |
| dsum | DSUM, PSUM, USUM | Pressure zone integrals (global) |
| edge | EDGE(:) | Edge viscous stress/vorticity |

## Read/Write Summary

| Routine | Reads | Writes |
|---------|-------|--------|
| COMPUTE_VISCOSITY | vel/vel_s, rho/rho_s, species/species_s, temp | visc (MU, KRES) |
| MASS_FINITE_DIFF | vel/vel_s, rho/rho_s, species/species_s, visc, temp | diff (DEL_RHO_D_DEL_Z) |
| DENSITY | vel/vel_s, rho/rho_s, species/species_s, diff, wall, pbar | rho_s/rho, species_s/species, temp |
| VELOCITY_FLUX | vel/vel_s, rho/rho_s, visc (MU), div/div_s, pres, edge | flux (FVX, FVY, FVZ + baroclinic) |
| WALL_BC | vel, rho/rho_s, species/species_s, temp, wall | temp(bdy), wall, pbar/pbar_s, energy (Q, M_DOT_PPP) |
| PARTICLE_MOMENTUM | part, vel/vel_s, flux (FVX) | flux (FVX += drag) |
| DIVERGENCE_PART_1 | vel/vel_s, rho/rho_s, species/species_s, temp, visc, energy, qr, pbar/pbar_s, diff, wall | div/div_s, dsum |
| DIVERGENCE_PART_2 | div/div_s, pbar/pbar_s, dsum, vel, vel_s, wall | div/div_s (final), dddt, pbar/pbar_s (D_PBAR_DT) |
| PRESSURE_SOLVE | flux, dddt, vel/vel_s, rho/rho_s, wall, visc (KRES) | pres/pres_s (H/HS) |
| VEL_PREDICTOR | vel, pres (H), flux | vel_s (US, VS, WS) |
| VEL_CORRECTOR | vel, vel_s, pres_s (HS), flux | vel (U, V, W) |
| COMBUSTION | species, temp, rho | species (ZZ), energy (Q, M_DOT_PPP) |
| PART_MASS_ENERGY | part, rho, species, temp | energy (M_DOT_PPP, Q), part |
| RADIATION | temp, rho, species, wall | qr (QR) |

## Pipelining Opportunities

### Predictor: Two-Level Pipeline

After MESH_EXCHANGE(1), the predictor currently runs VELOCITY_FLUX -> WALL_BC -> PARTICLE_MOMENTUM -> DIV_PART_1 sequentially. Data-flow analysis reveals two levels of independent execution:

**Level 1: VELOCITY_FLUX || WALL_BC**

These routines share no write conflicts:

| Array | VELOCITY_FLUX | WALL_BC | Conflict? |
|-------|---------------|---------|-----------|
| vel (U,V,W) | Read | Read | None (R\|\|R) |
| rho/rho_s | Read | Read | None (R\|\|R) |
| visc (MU) | Read | Read | None (R\|\|R) |
| flux (FVX,FVY,FVZ) | **Write** | -- | None |
| species/species_s | -- | Read | None |
| temp (TMP) | -- | Read/Write(bdy) | None |
| energy (Q, M_DOT_PPP) | -- | **Write** | None |
| pbar (D_PBAR_DT_S) | -- | **Write** | None |

VELOCITY_FLUX writes only to flux arrays. WALL_BC writes only to temp (boundary cells), wall, pbar, and energy. Zero overlap.

**Level 2: PARTICLE_MOMENTUM || DIVERGENCE_PART_1**

After Level 1 completes:
- PARTICLE_MOMENTUM reads FVX (from VELOCITY_FLUX) and adds drag forces (FVX += drag)
- DIVERGENCE_PART_1 reads energy and pbar (from WALL_BC) but does **NOT** read FVX/FVY/FVZ

These are independent because:
- PART_MOM touches flux (FVX) -- DIV_P1 does not
- DIV_P1 touches energy, pbar, div_s, dsum -- PART_MOM does not

**Join point**: PRESSURE_SOLVE, which needs both FVX (with drag, from PART_MOM) and DDDT (from DIV_P2, which follows DIV_P1).

### Corrector: Major Pipelining Opportunity

After MESH_EXCHANGE(4), the corrector has the biggest pipelining opportunity:

**Branch A (short)**: VELOCITY_FLUX -- writes FVX/FVY/FVZ, completes quickly.

**Branch B (long)**: COMBUSTION -> CONDENSATION -> PART_MASS_ENERGY -> MOVE_PARTICLES -- a chain of 4+ operations that takes much longer.

These are independent because:
- VELOCITY_FLUX writes flux; COMBUSTION/CONDENSATION/PME/MOVE don't touch flux
- COMBUSTION writes species/energy; VELOCITY_FLUX doesn't touch species/energy

**Join point**: PARTICLE_MOMENTUM, which needs both:
- FVX from Branch A (VELOCITY_FLUX)
- FVX_D and particle positions from Branch B (MOVE_PARTICLES)

After the join, execution continues sequentially: EX(7) -> WALL_BC -> EX(6) -> RADIATION -> EX(2) -> DIV_P1 -> DIV_EX -> DIV_P2 -> PRESSURE_SOLVE -> VEL_CORRECTOR.

The corrector opportunity is the bigger win: Branch B is a long chain, so VELOCITY_FLUX completes and its result sits ready well before PARTICLE_MOMENTUM needs it. This is essentially "free" overlap.

## Current Pipeline Execution Order

### Predictor Phase

```
INSERT_PARTICLES
  -> COMPUTE_VISCOSITY -> MASS_FINITE_DIFF -> DENSITY
  -> MESH_EXCHANGE(1)
  -> [pipeline start]
     Branch A: VELOCITY_FLUX -> PARTICLE_MOMENTUM ---------> PRESSURE_SOLVE
     Branch B: WALL_BC -> DIVERGENCE_PART_1 -> DIV_EX -> DIV_P2 -/
  -> VEL_PREDICTOR -> CHECK_STABILITY -> MESH_EXCHANGE(3) -> VEL_BC
```

### Corrector Phase

```
COMPUTE_VISCOSITY -> MASS_FINITE_DIFF -> DENSITY
  -> MESH_EXCHANGE(4)
  -> [pipeline start]
     Branch A: VELOCITY_FLUX --------------------------------> PARTICLE_MOMENTUM
     Branch B: COMBUSTION -> CONDENSATION -> PME -> MOVE ---/
  -> EX(7) -> WALL_BC -> EX(6) -> RADIATION -> EX(2)
  -> DIV_P1 -> DIV_EX -> DIV_P2
  -> PRESSURE_SOLVE -> VEL_CORRECTOR -> EX(6) -> VEL_BC -> DUMP
```

## Hedgehog Implementation Strategy

Pipelining in Hedgehog requires a different approach from the current scatter-gather pattern. Instead of dispatching N mesh tokens through a single task, we need to **fork a single mesh token into two concurrent branches** and join them at a synchronization point.

### Approach: Type-Based Branch Routing

Hedgehog routes data by C++ type. To fork, the orchestrator emits two different types:
- `VelocityFluxWork` -> routed to VELOCITY_FLUX task
- `WallBCWork` -> routed to WALL_BC task

A collector state waits for both results before emitting downstream.

### Key Considerations

1. **Per-mesh fork-join**: Each mesh independently forks into branches. With N meshes and 2 branches, there are 2N concurrent work items.
2. **Branch completion tracking**: The join-point collector must match results from the same mesh (use mesh index).
3. **Thread budget**: Pipelining adds concurrency orthogonal to mesh parallelism. With N meshes and 2 branches, the thread pool must accommodate up to 2N simultaneous tasks.
4. **Barrier semantics**: Mesh exchanges remain barriers across all meshes. Pipelining only applies *between* exchanges.

## Array Lifetime and Scratch Analysis

A second pass of analysis classifies each data array by its **lifetime** (persistent vs ephemeral) and **visibility** (output-visible vs purely internal). This matters for pipelining because ephemeral, internal arrays can be duplicated per-branch to eliminate false sharing.

### Array Persistence Classification

**Persistent arrays** carry meaningful state across timesteps. They are the "real" simulation state:

| Array | Output-Visible | Notes |
|-------|----------------|-------|
| U, V, W (velocity) | Yes (U/V/W-VELOCITY) | Updated by VEL_CORRECTOR |
| RHO (density) | Yes (DENSITY) | Updated by DENSITY each corrector |
| ZZ (species) | Yes (mass fractions) | Updated by DENSITY each corrector |
| TMP (temperature) | Yes (TEMPERATURE) | Recomputed from equation of state |
| H, HS (pressure head) | Yes (H, HS) | Updated by PRESSURE_SOLVE |
| PBAR, D_PBAR_DT | Yes (PRESSURE) | Background pressure, persistent |
| DEL_RHO_D_DEL_Z (diff) | Restart only | Mixed: zeroed in DIV_P1 but old value saved in DENSITY |
| EDGE | Indirect (via velocities) | Set at initialization, persists |
| WALL | Yes (BNDF) | Boundary cell data |
| LAGRANGIAN_PARTICLE | Yes (PART files) | Particle state |

**Ephemeral arrays** are fully recomputed each predictor or corrector phase. Their previous values are never carried forward:

| Array | Zeroed/Overwritten By | Output-Visible | Lifespan |
|-------|----------------------|----------------|----------|
| FVX, FVY, FVZ (flux) | Overwritten by VELOCITY_FLUX | Yes (F_X/F_Y/F_Z) | VFLUX -> PRESSURE_SOLVE + VEL_PRED/CORR |
| FVX_D, FVY_D, FVZ_D (drag) | Explicitly zeroed by MOVE_PARTICLES | Yes (DRAG FORCE) | MOVE -> PART_MOM (very short) |
| D / DS (divergence) | Zeroed by DIV_P1 | Yes (DIVERGENCE) | DIV_P1 -> PRESSURE_SOLVE |
| **DDDT** | Recomputed by DIV_P2 | **No** | DIV_P2 -> PRESSURE_SOLVE only |
| MU (viscosity) | Overwritten by COMPUTE_VISCOSITY | Yes (VISCOSITY) | VISC -> VFLUX + DIV_P1 |
| KRES (kinetic energy) | Overwritten by COMPUTE_VISCOSITY | Yes (RESOLVED KE) | VISC -> PRESSURE_SOLVE |
| US, VS, WS (vel_s) | Overwritten by VEL_PREDICTOR | Partial (CFL) | VEL_PRED -> corrector input |
| RHOS (density_s) | Overwritten by DENSITY | No | DENSITY -> corrector input |
| ZZS (species_s) | Overwritten by DENSITY | No | DENSITY -> corrector input |
| Q (energy source) | Reset by COMBUSTION/WALL_BC | Yes (HRRPUV) | Written -> consumed by DIV_P1 |
| M_DOT_PPP | Zeroed in DENSITY | Restart only | PART_MASS_ENERGY -> DENSITY/DIV_P1 |
| QR (radiation) | Zeroed by RADIATION | Yes (RAD LOSS) | RADIATION -> DIV_P1 |
| **DSUM, PSUM, USUM** | Zeroed by InitDiv | **No** | DIV_P1 -> DIV_EXCHANGE -> DIV_P2 |
| **D_PBAR_DT_S** | Recomputed by DIV_P2 | **No** | DIV_P2 -> PRESSURE_SOLVE |
| PBAR_S | Computed from PBAR in DENSITY | No | DENSITY -> DIV_P1/P2 |

### Scratch Arrays (WORK)

MESH_TYPE contains shared scratch arrays reused across routines:

| Scratch Pool | Dimensions | Used By |
|-------------|------------|---------|
| WORK1-9 | 3D (0:IBP1, 0:JBP1, 0:KBP1) | VELOCITY_FLUX (1-6), COMPUTE_VISCOSITY (1-6), DIV_P1 (1-7,9), DENSITY (4-5), RADIATION (1-9), PART_MASS_ENERGY (1-2,4-7) |
| SWORK1-4 | 4D (+ N_SCALARS) | DIV_P1 (1-3), DENSITY (4), PART_MASS_ENERGY (1) |
| TURB_WORK1-10 | 3D | COMPUTE_VISCOSITY only |
| PRHS, BX\*/BY\*/BZ\* | Solver-specific | PRESSURE_SOLVE only |
| IWORK1 | 3D integer | COMPUTE_VISCOSITY only |
| WALL_WORK1-2 | 1D (N_WALL_CELLS) | RADIATION only |
| FACE_WORK1-3 | 1D | RADIATION only |

These arrays are not persistent -- they are overwritten at the start of each routine that uses them. They exist solely to avoid repeated allocation.

### WORK Array Conflict Matrix for Pipeline Candidates

Because WORK1-9 are shared per-mesh, two routines that use overlapping WORK arrays cannot run concurrently on the same mesh without corruption. This is a **hidden dependency** not visible in the logical data-flow graph:

| Parallel Pair | WORK Conflicts | Safe? |
|---------------|----------------|-------|
| VELOCITY_FLUX \|\| WALL_BC | None (WALL_BC uses no WORK arrays) | **Yes** |
| PARTICLE_MOMENTUM \|\| DIV_P1 | **Needs investigation** | **Unknown** |
| VELOCITY_FLUX \|\| COMBUSTION | None (COMBUSTION uses no WORK arrays) | **Yes** |
| VELOCITY_FLUX \|\| CONDENSATION | **Needs investigation** | **Unknown** |
| VELOCITY_FLUX \|\| PART_MASS_ENERGY | WORK1-2, 4-7 conflict | No (without mitigation) |
| VELOCITY_FLUX \|\| RADIATION | WORK1-9 all conflict | No (without mitigation) |
| DIV_P1 \|\| RADIATION | WORK1-9 all conflict | No (without mitigation) |

For our identified pipeline candidates:
- **Predictor Level 1** (VELOCITY_FLUX \|\| WALL_BC): **Confirmed safe** -- no WORK conflicts.
- **Predictor Level 2** (PARTICLE_MOMENTUM \|\| DIV_P1): PARTICLE_MOMENTUM_TRANSFER in part.f90 does not use WORK arrays directly, but **further investigation is needed** to confirm no indirect WORK usage through called subroutines.
- **Corrector** (VELOCITY_FLUX \|\| COMBUSTION->COND->PME->MOVE): COMBUSTION and CONDENSATION don't use WORK arrays, so they are safe to overlap with VELOCITY_FLUX. PART_MASS_ENERGY uses WORK1-2,4-7, but VELOCITY_FLUX is expected to complete before PME begins (it is on the short branch). However, **if VELOCITY_FLUX is slow for a particular mesh, a timing-dependent race is possible** -- this needs a synchronization guarantee or mitigation.

### Mitigation: Per-Branch Scratch Allocation

The WORK array conflicts can be eliminated entirely by allocating **separate scratch pools per pipeline branch**. Since these arrays are:
- Ephemeral (overwritten at the start of each routine)
- Not output-visible (never appear in result files)
- Not persistent (no state carries between routines)

They can be duplicated without affecting simulation correctness.

**Approach**: Allocate a second set of scratch arrays (WORK1B-9B, SWORK1B-4B, etc.) on each mesh. Each pipeline branch uses its own pool:

```
Branch A (VELOCITY_FLUX path):  uses WORK1-9 (original)
Branch B (WALL_BC/DIV_P1 path): uses WORK1B-9B (new)
```

This eliminates all WORK conflicts and opens up additional pipelining candidates that were blocked by false sharing. The memory cost is modest -- 9 3D arrays per mesh (~72 bytes/cell for double precision, or ~4.5 MB for a 64^3 mesh).

The same principle applies to any purely internal ephemeral array that creates a false dependency between branches. Candidates for duplication:

| Array | Duplicable? | Reason |
|-------|-------------|--------|
| WORK1-9 | Yes | Pure scratch, no output, no persistence |
| SWORK1-4 | Yes | Pure scratch for species computations |
| TURB_WORK1-10 | Yes | Pure scratch for viscosity model |
| PRHS | Yes | Pressure solver internal |
| BXS/BXF/BYS/BYF/BZS/BZF | Yes | Pressure solver boundary conditions |
| DDDT | Yes | Never output, only flows DIV_P2 -> PRESSURE_SOLVE |
| D_PBAR_DT_S | Yes | Never output, internal pressure derivative |
| DSUM/PSUM/USUM (local) | Yes | Per-mesh accumulators, zeroed each phase |

**Not duplicable** (true shared state read by multiple branches):

| Array | Why Not |
|-------|---------|
| U, V, W | Read by both VELOCITY_FLUX and WALL_BC |
| RHO, RHOS | Read by both branches |
| TMP | Read by both branches |
| ZZ, ZZS | Read by both branches |
| MU | Read by both VELOCITY_FLUX and DIV_P1 |
| PBAR_S | Read by both WALL_BC and DIV_P1 |

These shared reads are safe (R\|\|R) and don't need duplication.

### Expanded Pipeline Candidates (With Scratch Duplication)

With per-branch scratch pools, the following additional pipeline candidates become feasible:

| Parallel Pair | Previously Blocked By | Status |
|---------------|----------------------|--------|
| DIV_P1 \|\| RADIATION | WORK1-9 conflict | **Feasible with scratch duplication** -- but logical data dependency needs verification (does DIV_P1 read QR from RADIATION?) |
| VELOCITY_FLUX \|\| PART_MASS_ENERGY | WORK1-2,4-7 conflict | **Feasible with scratch duplication** -- but only useful if VELOCITY_FLUX hasn't already completed |
| VELOCITY_FLUX \|\| RADIATION | WORK1-9 conflict | **Feasible with scratch duplication** -- but logical data flow must be verified |

**Important**: Scratch duplication removes the *mechanical* conflict but does not override *logical* data dependencies. Each candidate above still requires verification that the routines don't share logical read-write dependencies on physics arrays (the Read/Write Summary table above). Scratch duplication only helps when the sole blocking dependency was the shared WORK arrays.

### Items Requiring Further Investigation

1. **PARTICLE_MOMENTUM WORK usage**: Verify that PARTICLE_MOMENTUM_TRANSFER and its callees in part.f90 do not use any WORK arrays through indirect calls. If confirmed clean, predictor Level 2 parallelism is fully safe.

2. **CONDENSATION scratch usage**: Determine whether CONDENSATION (fire.f90) uses WORK arrays. If not, the full corrector Branch B chain (COMBUSTION->CONDENSATION->PME->MOVE) is clean for overlap with VELOCITY_FLUX without needing scratch duplication.

3. **Corrector timing guarantee**: In the corrector pipeline, VELOCITY_FLUX (Branch A) is expected to finish before PART_MASS_ENERGY (Branch B) starts using WORK arrays. This assumption depends on VELOCITY_FLUX being faster than COMBUSTION+CONDENSATION combined. If this timing assumption is unreliable, scratch duplication is needed as a safety measure.

4. **EDGE array thread safety**: EDGE is persistent and read by VELOCITY_FLUX. If any routine on the parallel branch writes to EDGE, this creates a hidden conflict. Needs verification that WALL_BC and DIV_P1 do not modify EDGE.

5. **DEL_RHO_D_DEL_Z history dependency**: This array is zeroed in DIV_P1 but its previous value is saved by DENSITY (into SWORK4) at the start of the phase. Since the save happens before the pipeline fork, both branches can safely read the saved copy. Needs confirmation that no branch writes to DEL_RHO_D_DEL_Z before DIV_P1.

## Source Files Analyzed

- `Source/main.f90` -- predictor/corrector phase sequencing
- `Source/velo_kernels.f90` -- VELOCITY_FLUX_KERNEL, VEL_PREDICTOR/CORRECTOR_KERNEL
- `Source/divg_kernels.f90` -- DIVERGENCE_PART_1_KERNEL, DIVERGENCE_PART_2_KERNEL
- `Source/pres_kernels.f90` -- PRESSURE_SOLVER_COMPUTE_RHS, FFT solve
- `Source/wall_kernels.f90` -- WALL_BC kernels
- `Source/part.f90` -- PARTICLE_MOMENTUM_TRANSFER, MOVE_PARTICLES
- `Source/fire.f90` -- COMBUSTION
- `Source/radi.f90` -- COMPUTE_RADIATION
- `Source/turb.f90` -- COMPUTE_VISCOSITY
- `Source/mass.f90` -- MASS_FINITE_DIFFERENCES, DENSITY
