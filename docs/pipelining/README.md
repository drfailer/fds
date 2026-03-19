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
