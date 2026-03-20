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
| [fds_section_pipeline.dot](fds_section_pipeline.dot) / [.svg](fds_section_pipeline.svg) | Section-level pipeline structure with cost annotations and Hedgehog graph recommendations |

**Implementation**: [IMPLEMENTATION_PROGRESS.md](IMPLEMENTATION_PROGRESS.md) — phased implementation plan and progress tracking.

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

### Investigation Results

Items 1-2 from the previous section have been resolved by the intra-routine analysis below:

1. **PARTICLE_MOMENTUM WORK usage**: **Resolved -- CONFIRMED SAFE.** PARTICLE_MOMENTUM_TRANSFER_KERNEL (part_kernels.f90) uses NO WORK arrays. It only accesses FVX/FVY/FVZ, FVX_D/FVY_D/FVZ_D, and velocity arrays. Predictor Level 2 parallelism (PART_MOM || DIV_P1) is safe from scratch conflicts.

2. **CONDENSATION scratch usage**: **Resolved -- USES WORK1-2 and SWORK1.** CONDENSATION_EVAPORATION_KERNEL uses WORK1 (RHO_INTERIM), WORK2 (TMP_INTERIM), and SWORK1 (ZZ_INTERIM) as snapshots. However, the per-cell computation is independent and these are read-only snapshots taken at the start. If VELOCITY_FLUX completes before CONDENSATION begins (expected since COMBUSTION runs first on Branch B), there is no conflict. If timing cannot be guaranteed, scratch duplication is needed.

3. **Corrector timing guarantee**: Still needs runtime profiling. COMBUSTION uses no WORK arrays, so VELOCITY_FLUX can safely overlap with COMBUSTION. CONDENSATION starts after COMBUSTION and uses WORK1-2/SWORK1, creating a potential conflict if VELOCITY_FLUX hasn't completed. For safety, scratch duplication is recommended.

4. **EDGE array thread safety**: EDGE is read-only during the pipeline window. VELOCITY_FLUX reads EDGE for vorticity/stress interpolation. Neither WALL_BC nor DIV_P1 write to EDGE. EDGE is only modified at initialization (init.f90) and during obstruction creation/removal (REDEFINE_EDGE). **Confirmed safe.**

5. **DEL_RHO_D_DEL_Z history dependency**: The old value is saved into SWORK4 by DENSITY before the pipeline fork. DIV_P1 zeroes and recomputes DEL_RHO_D_DEL_Z during its execution. No other branch writes to it before DIV_P1. **Confirmed safe.**

## Intra-Routine Decomposition Analysis

Beyond pipelining entire routines, we analyzed the internal structure of each expensive routine to identify independent sections that could be split for finer-grained parallelism.

### VELOCITY_FLUX -- Two-Phase Structure

VELOCITY_FLUX_KERNEL (velo_kernels.f90:328-916) has a clear two-phase structure:

**Phase 1: Shared intermediates** (~35 lines)
- Computes vorticity (OMX, OMY, OMZ) and viscous stress tensor (TXY, TXZ, TYZ) over the full grid
- Uses WORK1-6 as scratch for these 6 fields
- Must complete before Phase 2

**Phase 2: Three independent flux computations**

| Section | Loop Range | Reads | Writes | Independent? |
|---------|-----------|-------|--------|--------------|
| FVX | K=1:KBAR, J=1:JBAR, I=0:IBAR | OMY, OMZ, TXY, TXZ, MU, VV, WW, EDGE | M%FVX | Yes (after Phase 1) |
| FVY | K=1:KBAR, J=0:JBAR, I=1:IBAR | OMX, OMZ, TXY, TYZ, MU, UU, WW, EDGE | M%FVY | Yes (after Phase 1) |
| FVZ | K=0:KBAR, J=1:JBAR, I=1:IBAR | OMX, OMY, TXZ, TYZ, MU, UU, VV, EDGE | M%FVZ | Yes (after Phase 1) |

FVX, FVY, FVZ write to distinct arrays and read shared intermediates (read-only after Phase 1). All EDGE accesses are reads. **These three can run in parallel.**

**BAROCLINIC_CORRECTION_KERNEL** (velo_kernels.f90:31-108): After Phase 1 setup (P, RRHO using WORK1-2), the FVX_B, FVY_B, FVZ_B corrections are fully independent of each other.

**Decomposition opportunity**: Split VELOCITY_FLUX into:
1. Vorticity + stress tensor (shared, sequential)
2. FVX, FVY, FVZ (independent, parallel)

**Estimated benefit**: Modest. The three loops have similar cost and the shared Phase 1 is ~30% of total. With 3-way split: theoretical 1.5x within VELOCITY_FLUX.

### COMPUTE_VISCOSITY -- Independent MU and KRES

COMPUTE_VISCOSITY_KERNEL (velo_kernels.f90:1270-1706) has two independent outputs:

| Section | Reads | Writes | Dependencies |
|---------|-------|--------|-------------|
| MU_DNS | TMP, ZZ | M%MU_DNS | None |
| STRAIN_RATE | Velocities | M%STRAIN_RATE | None (wall loop sequential) |
| Turbulent MU | MU_DNS, STRAIN_RATE, RHO | M%MU | After MU_DNS + STRAIN_RATE |
| **KRES** | **UU, VV, WW only** | **M%KRES** | **None -- fully independent** |
| Wall mirroring | MU, KRES | M%MU, M%KRES (ghost) | After both MU and KRES |

**KRES can run in parallel with the entire MU computation chain.** KRES reads only velocity arrays and writes only to M%KRES. No dependency on MU_DNS, STRAIN_RATE, or the turbulence model.

**Limitation**: DEARDORFF and DYNSMAG turbulence models require FILL_EDGES and TEST_FILTER (global operations) and cannot be K-decomposed. CONSMAG, VREMAN, WALE are fully K-decomposable.

### DIVERGENCE_PART_1 -- Six Phases with Barriers

DIVERGENCE_PART_1_KERNEL (divg_kernels.f90:28-1396) is the largest and most complex routine. It has 6 major phases with internal dependencies:

```
Phase 1: Setup (zero DP, pointer aliases)
  |
  v
Phase 2A: Species diffusion fluxes (RHO_D_DZDX/Y/Z)
  |  Mass conservation barrier (MAXLOC/SUM across species)
  v
Phase 2B: Diffusive heat flux (H_RHO_D_DZDX/Y/Z -> DP)
  |  Depends on Phase 2A output
  |
  +--- Phase 3: Specific heat (CP, R_H_G)  [INDEPENDENT of Phase 2]
  |      |
  |      v
  |    Phase 4: Thermal conductivity (KP) -> Thermal divergence (-> DP)
  |      |  Depends on Phase 3 (conditional)
  |
  +--- Phase 5A: Enthalpy advection (-> DP)  [INDEPENDENT of Phases 2-4]
  |
  v  (all phases accumulate into DP)
Phase 5B: RTRM = 1/(rho*CP*TMP)  [depends on Phase 3 if !CONSTANT_SPECIFIC_HEAT_RATIO]
  |  DP *= RTRM (multiplicative scaling of all accumulated terms)
  v
Phase 5C: Species advection Part 1 (FX_ZZ, FY_ZZ, FZ_ZZ face fluxes)
  |  MW correction barrier (MAXLOC/SUM across species)
  v
Phase 5D: Species advection Part 2 (per-species divergence -> DP)
  |  Depends on Phase 5C output
  v
Phase 5E-G: Source terms (reactions, stratification, MMS -> DP)
  v
Phase 6: Pressure zone sums (DSUM, PSUM, USUM from final DP)
```

**Independent sections that can run in parallel:**

| Parallel Group | Sections | Constraint |
|---------------|----------|------------|
| Group A | Phase 2A+2B (species/heat diffusion) | Sequential internally (2A -> 2B) |
| Group B | Phase 3+4 (specific heat + thermal conductivity + thermal divergence) | Sequential internally (3 -> 4) |
| Group C | Phase 5A (enthalpy advection) | Independent of A and B |
| **A \|\| B \|\| C** | All three groups | **Yes, can run in parallel** |

After groups A, B, C complete and their contributions are accumulated into DP, Phase 5B applies the RTRM scaling, then Phases 5C-6 run sequentially.

**Barriers preventing further decomposition:**
- Phase 2A: Mass conservation correction requires all species fluxes (global reduction per cell)
- Phase 5C: MW correction requires all species face fluxes (global reduction per face)
- Phase 5B: DP *= RTRM is a full-array multiplicative gate between additive accumulation (2B+4+5A) and species advection (5C+5D)
- Phase 6: Global accumulation into DSUM/PSUM/USUM

**Estimated benefit**: Groups A, B, C are roughly equal cost (~200 lines each). Running them in parallel could yield ~2-3x speedup within DIV_P1. However, this requires separate scratch arrays for each group (they all use WORK arrays differently).

### RADIATION -- Angle-Level Parallelism (Major Opportunity)

COMPUTE_RADIATION (radi.f90, ~1200 lines) has the biggest intra-routine parallelism opportunity:

**Phase 1: Absorption coefficients** (~350 lines)
- Computes KAPPA_GAS, KFST4_GAS, EXTCOE, KAPPA_PART, SCAEFF per cell
- K-block decomposable (per-cell independent)
- Must complete before angle loop

**Phase 2: Angle loop** (~430 lines) -- **MAJOR OPPORTUNITY**
- Sweeps NUMBER_RADIATION_ANGLES angles (typically 100-500)
- Only ANGLE_INCREMENT angles updated per radiation call
- Each angle N:
  1. Set boundary intensity IL from wall data for angle N (independent)
  2. Sweep cells in upwind order for angle N (sequential within angle)
  3. Update UII accumulator (reduction across angles)
  4. Update wall outgoing intensity (independent per angle)

**Key insight**: Each angle sweep is **fully independent** of other angles. The cell sweep within each angle has sequential I->J->K upwind dependency, but different angles can run simultaneously.

| Aspect | Detail |
|--------|--------|
| Parallelism grain | NUMBER_RADIATION_ANGLES (100-500) |
| Per-angle private data | IL (WORK2), IL_UP (WORK8) |
| Shared read-only data | KFST4_GAS (WORK1), EXTCOE (WORK4), KAPPA_PART (WORK5), SCAEFF (WORK6), KFST4_PART (WORK7) |
| Reduction at end | UIID accumulator (+=), INRAD_W wall flux (+=) |
| Memory cost per thread | 2 3D arrays (IL, IL_UP) |

**Estimated benefit**: With 10 angle threads on a 100-angle problem, ~10x speedup for the angle loop (the dominant cost of RADIATION). This is the single largest parallelism opportunity in FDS.

**Phase 3: QR assembly** (~50 lines)
- QR = KAPPA_GAS * UII - KFST4_GAS (per-cell, K-decomposable)
- Must run after all angle sweeps complete

### COMBUSTION -- Cell-Level Parallelism

COMBUSTION_GENERAL_KERNEL (fire.f90:532+) has two phases:

**Phase 1: Identify active cells** (~50 lines)
- Filters cells by species/temperature thresholds
- Builds sparse active cell list
- K-decomposable

**Phase 2: Chemistry ODE per cell** (~30 lines of loop, but COMBUSTION_MODEL is expensive)
- Each cell solves an independent ODE system (species + energy)
- **Embarrassingly parallel** across active cells
- No cell-to-cell coupling
- Uses NO WORK arrays

Reactions within a single cell are NOT parallelizable (coupled ODE system solved by CVODE or fast chemistry). But the cell-level parallelism is the right grain -- active cells are typically 0-10% of the mesh, and each COMBUSTION_MODEL call is expensive (100-1000+ FLOPs).

### DENSITY -- Species-Level Parallelism

DENSITY_KERNEL (mass_kernels.f90:351-920) has one key parallelizable section:

**Species advection** (N=1:N_TOTAL_SCALARS outer loop):
- Each species N computes M%ZZS(:,:,:,N) independently
- No cross-species dependency within the advection loop
- Can split species across threads

**Barriers**:
- Density summation: RHOS = SUM(ZZS, dim=species) -- requires all species complete
- CHECK_MASS_DENSITY: 6-neighbor scatter -- forces sequential post-processing
- Temperature update: TMP = PBAR / (RSUM * RHOS) -- after density summation

**Estimated benefit**: With N_TOTAL_SCALARS species (typically 3-10), species-level parallelism gives 3-10x for the advection loop. Already handled by the existing DENSITY_BLOCK_KERNEL K-decomposition.

### WALL_BC -- Per-Wall-Cell and Per-Species Parallelism

WALL_BC kernels (wall_kernels.f90) are organized as independent per-wall-cell operations:

| Sub-kernel | Parallelism | Notes |
|-----------|-------------|-------|
| NEAR_SURFACE_GAS_VARIABLES | Per wall cell | Each wall cell reads gas-phase data independently |
| CALCULATE_RHO_F | Per wall cell | Each wall cell computes surface density independently |
| ASSIGN_GHOST_VALUE | Per wall cell | Each external wall writes to distinct ghost cell |
| CALC_DEPOSITION | Per species | N_TRACKED_SPECIES independent deposition velocity calculations |
| PYROLYSIS | Per material | N_MATS independent reaction rates, then sequential accumulation |
| SOLID_HEAT_TRANSFER | Per wall cell | 1D conduction solve per wall cell (expensive, independent) |

**SOLID_HEAT_TRANSFER is the expensive part** -- it solves a 1D heat equation through the solid for each wall cell. These are fully independent and embarrassingly parallel.

Uses NO WORK arrays (confirmed earlier). Per-species and per-material loops offer secondary parallelism.

### PRESSURE_SOLVE -- Sequential FFT Core

PRESSURE_SOLVER_COMPUTE_RHS + FFT (pres_kernels.f90):

| Section | Parallelizable? | Notes |
|---------|----------------|-------|
| Boundary condition setup (BXS/BXF/BYS/BYF/BZS/BZF) | Yes, per wall cell | Non-overlapping boundary arrays by IOR direction |
| PRHS computation | Yes, per (I,J,K) | Pure divergence of flux terms |
| **FFT solve** | **No** | **Inherently sequential** (global transform) |
| Solution copy to H/HS | Yes, per (I,J,K) | Simple array copy |
| H boundary conditions | Yes, per face | Independent per boundary face |

The FFT solve is the bottleneck. ULMAT (sparse direct solver) is also inherently sequential per pressure zone but can parallelize across zones. The RHS and boundary setup (~30% of total) can be parallelized.

### PARTICLE_MASS_ENERGY -- Particle Loop (Not Decomposable)

PARTICLE_MASS_ENERGY_KERNEL (part.f90:3464-4588):

- Main PARTICLE_LOOP iterates over all particles sequentially
- Each particle's heat/mass transfer is independent of other particles
- **But**: particles accumulate into shared grid arrays (M_DOT_PPP, Q at cell I,J,K)
- Multiple particles in the same cell create write conflicts
- Uses WORK1-7 and SWORK1

**Parallelism**: Could use atomic accumulation or particle-cell binning, but the current structure mixes particle classes in one loop. Not easily decomposable without restructuring.

### MOVE_PARTICLES -- Sequential (Mesh Transfer)

MOVE_PARTICLES (part.f90:1807-3456):

- Uses NO WORK arrays
- Sequential due to inter-mesh particle transfers and shared LAGRANGIAN_PARTICLE array
- Particle removal invalidates indices, preventing parallel iteration
- Not a candidate for intra-routine parallelism

### Summary: Intra-Routine Parallelism Opportunities

Ranked by estimated impact:

| Rank | Routine | Opportunity | Parallelism Type | Est. Speedup | Complexity |
|------|---------|-------------|-----------------|--------------|------------|
| 1 | **RADIATION** | Angle loop (100-500 angles) | Angle-parallel | **10-50x** for angle loop | Medium (private IL per thread, UIID reduction) |
| 2 | **DIV_P1** | Groups A\|\|B\|\|C (diffusion \|\| thermal \|\| enthalpy advection) | Section-parallel | **2-3x** within DIV_P1 | High (separate WORK pools, DP accumulation sync) |
| 3 | **COMBUSTION** | Per-cell chemistry ODE | Cell-parallel | **Nx** (N = active cells / threads) | Low (embarrassingly parallel, no WORK) |
| 4 | **VELOCITY_FLUX** | FVX \|\| FVY \|\| FVZ after vorticity | Component-parallel | **~1.5x** within VFLUX | Low (distinct output arrays) |
| 5 | **DENSITY** | Per-species advection | Species-parallel | **3-10x** for advection loop | Low (independent 4D slices) |
| 6 | **WALL_BC** | SOLID_HEAT_TRANSFER per wall cell | Wall-cell-parallel | **Nx** (N = wall cells / threads) | Low (independent 1D solves) |
| 7 | **COMPUTE_VISCOSITY** | KRES \|\| MU chain | Section-parallel | **~1.3x** within VISC | Low (independent outputs) |
| 8 | **PRESSURE_SOLVE** | RHS + boundary setup | Cell-parallel | **~1.3x** (30% of solve) | Low (FFT is bottleneck) |

**Not decomposable**: PARTICLE_MASS_ENERGY (particle-cell accumulation conflicts), MOVE_PARTICLES (mesh transfers, index invalidation), FFT solve (global transform).

## Operation Reordering Analysis

By examining data dependencies at the section level (not just routine level), we can identify reordering opportunities that extend the parallel windows beyond what routine-level pipelining achieves.

### Key Discovery: Interior vs Wall Corrections

Most DIV_P1 sections have a two-part structure:
1. **Bulk interior loops** (triple I,J,K) -- read only mesh arrays (TMP, RHO, ZZ, MU) set by DENSITY/COMPUTE_VISCOSITY before the pipeline fork
2. **Wall correction loops** (iterate over wall cells) -- read B1 boundary properties (TMP_F, RHO_F, ZZ_F, U_NORMAL_S, RHO_D_F) set by WALL_BC

The interior computation does NOT need WALL_BC output. Only the wall correction loops depend on WALL_BC.

> **INVESTIGATION RESULT (Phase 6):** While individual wall correction loops are self-contained, the interior/wall split is **NOT feasible** for pipelining. The interior I,J,K loops include boundary face positions (e.g., I=0..IBAR), and the divergence accumulation reads these boundary values into DP. Crucially, Phase 5's `DP *= RTRM` is a **multiplicative** operation — any error from uncorrected boundary face values is amplified and cannot be fixed by an additive correction afterward. A clean split would require re-running the entire divergence computation for wall-adjacent cells, negating any savings. See [IMPLEMENTATION_PROGRESS.md](IMPLEMENTATION_PROGRESS.md) Phase 6 for the full analysis.

### WALL_BC Exact Writes (Predictor Phase)

WALL_BC writes during the predictor:

| Target | Arrays | Notes |
|--------|--------|-------|
| Ghost cells only | M%TMP, M%RHOS, M%ZZS at (BC%II, BC%JJ, BC%KK) | Never interior cells |
| B1 structure | B1%TMP_F, B1%RHO_F, B1%ZZ_F, B1%RHO_D_F, B1%U_NORMAL_S, B1%Q_RAD_OUT, B1%Q_CON_F | Wall surface properties |
| **NOT written** | M%Q, M%M_DOT_PPP, M%PBAR_S, M%D_PBAR_DT_S | Q only by COMBUSTION (corrector); PBAR_S only by DENSITY |

### DIV_P1 Section Dependencies on WALL_BC

| DIV_P1 Section | Interior (bulk I,J,K) | Wall Corrections | Depends on WALL_BC? |
|----------------|----------------------|------------------|---------------------|
| Phase 2A: Species diffusion fluxes | Reads ZZP, RHOP, TMP (interior) | Reads B1%RHO_D_DZDN_F | Interior: **No**. Wall: **Yes** |
| Phase 2B: Diffusive heat flux | Reads TMP, RHO_D_DZDX/Y/Z | Reads B1%TMP_F, B1%U_NORMAL_S, B1%RHO_D_DZDN_F | Interior: **No**. Wall: **Yes** |
| Phase 3: Specific heat (CP) | Reads ZZP, TMP | No wall loop | **No** |
| Phase 4: Thermal conductivity + divergence | Reads TMP, MU, ZZP | Reads B1%K_G (conditional) | Interior: **No**. Wall: **Partial** |
| Phase 5A: Enthalpy advection | Reads RHOP, TMP, ZZP, velocities | Reads B1%U_NORMAL_S, B1%TMP_F, B1%ZZ_F, B1%RHO_F | Interior: **No**. Wall: **Yes** |
| Phase 5B: RTRM | Reads R_H_G, RHOP | No wall loop | **No** |
| Phase 5C: Species advection Part 1 | Reads RHOP, ZZP, velocities | Metadata only (IOR checks) | **No** |
| Phase 5D: Species advection Part 2 | Reads FX_ZZ/FY_ZZ/FZ_ZZ | Reads B1%U_NORMAL_S, B1%RHO_F, B1%ZZ_F | Interior: **No**. Wall: **Yes** |
| Phase 5E-G: Source terms | Reads M%D_SOURCE, velocities | No wall loop | **No** |
| Phase 6: Pressure zone sums | Reads DP, RTRM | Reads B1%U_NORMAL_S (USUM accumulation) | Interior: **No**. Wall: **Yes** |

### Reordering Opportunity 1: VFLUX || WALL_BC (Predictor)

> **Note:** A three-way fork (VFLUX || WALL_BC || DIV_P1 interior) was investigated and ruled out. See the investigation note above. The feasible option is a simpler two-way fork.

Run VELOCITY_FLUX and WALL_BC concurrently after MESH_EXCHANGE(1). No data conflicts (confirmed in dependency proof table). PARTICLE_MOMENTUM chains after VFLUX (needs FVX). DIV_P1 runs sequentially after the join.

```
MESH_EXCHANGE(1)
├── Branch A: VELOCITY_FLUX → PARTICLE_MOMENTUM (266 ops/cell)
└── Branch B: WALL_BC (65 ops/cell, hidden behind A)

→ Join → DIV_P1 (full, 1300 ops/cell) → DIV_EXCHANGE → ...
```

**Estimated benefit**: 65 ops/cell saved (WALL_BC hidden behind VFLUX). Modest but reliable.

**Implementation complexity**: Low. Reuses the fork/join pattern from the corrector. No Fortran kernel changes needed.

### Reordering Opportunity 2: Pressure BC Setup || DIV_P2 (Predictor & Corrector)

PRESSURE_SOLVER_COMPUTE_RHS has two independent sections:

| Section | Reads | Needs DDDT? |
|---------|-------|-------------|
| Wall loop: BXS/BXF/BYS/BYF/BZS/BZF setup (lines 57-226) | H/HS, FVX/FVY/FVZ (boundary only), WALL_WORK1, KRES, velocities | **No** |
| PRHS computation (lines 231-298) | FVX/FVY/FVZ (all cells), **DDDT** | **Yes** |

The boundary condition setup reads H/HS (from previous iteration), FVX at boundary faces, and wall data. **It does not read DDDT at all.** This means:

```
Current:  DIV_P2 ──→ PRESSURE_SOLVE (BC setup + PRHS + FFT)

Reordered:
├── DIV_P2 ────────→ PRHS computation ──→ FFT solve
└── Pressure BC setup (BXS/BXF/...)  ───/
    [runs in parallel with DIV_P2]
```

**Estimated benefit**: Modest. BC setup is ~20-30% of pre-solve work. But since DIV_P2 is relatively fast, the overlap window is small.

### Reordering Opportunity 3: DIV_P1 Non-QR Sections || RADIATION (Corrector)

In the corrector, RADIATION writes QR and DIV_P1 reads it. But QR is consumed in only **one line** of DIV_P1:

```fortran
! divg_kernels.f90, line 566/577 in COMPUTE_THERMAL_DIVERGENCE:
DP(I,J,K) = DP(I,J,K) + DELKDELT + M%Q(I,J,K) + M%QR(I,J,K)
```

All other DIV_P1 sections (species diffusion, enthalpy advection, species advection, pressure zone sums) do NOT read QR. This means ~85% of DIV_P1 can start before RADIATION completes:

```
Corrector current:
  ... → WALL_BC → EX(6) → RADIATION → EX(2) → DIV_P1 (full) → ...

Corrector reordered:
  ... → WALL_BC → EX(6)
  ├── RADIATION (computes QR)
  └── DIV_P1 non-QR sections (~85% of work)
      ← Both complete
      DIV_P1 QR addition (1 line: DP += Q + QR)
      → DIV_EXCHANGE → DIV_P2 → ...
```

**Estimated benefit**: Significant if RADIATION is expensive (it often is). The non-QR sections of DIV_P1 run "for free" behind RADIATION.

**Note**: EX(2) currently sits between RADIATION and DIV_P1. EX(2) exchanges QR across meshes. The non-QR sections of DIV_P1 don't need QR, so they can start before EX(2). Only the QR addition line needs to wait for EX(2).

### Reordering Opportunity 4: Extended Corrector Three-Way Fork

Combining the corrector pipeline with the DIV_P1 || RADIATION overlap:

```
MESH_EXCHANGE(4)
├── Branch A: VELOCITY_FLUX (short, writes FVX)
├── Branch B: COMBUSTION → CONDENSATION → PME → MOVE (long chain)
│
│   ← Both branches join at PARTICLE_MOMENTUM
│
│   → EX(7) → WALL_BC → EX(6)
│   ├── Branch C: RADIATION (writes QR)
│   └── Branch D: DIV_P1 non-QR sections (species diffusion, enthalpy, species advection, ...)
│
│   ← Both branches join
│
│   DIV_P1 QR addition → DIV_EXCHANGE → DIV_P2
│   ├── Branch E: PRHS computation (needs DDDT from DIV_P2)
│   └── Branch F: Pressure BC setup (BXS/BXF, independent of DDDT)
│
│   ← Both branches join
│
│   FFT solve → VEL_CORRECTOR → ...
```

This creates **three fork-join pairs** in the corrector, each overlapping independent computations.

### Rejected Reorderings

| Candidate | Why Rejected |
|-----------|-------------|
| WALL_BC during VELOCITY_FLUX Phase 1 | MESH_EXCHANGE(1) is a hard barrier between them. WALL_BC needs exchanged boundary data. |
| Start DIV_P1 Phase 3 (CP) as separate early task | Fused with PARTICLE_MOMENTUM in current kernel. Overhead of splitting exceeds gain. |
| VEL_PREDICTOR BC setup before PRESSURE_SOLVE | H gradient `(H(I+1)-H(I))` is intrinsic to the velocity update formula. Cannot separate. |
| Move VELOCITY_FLUX earlier in corrector | Already optimally positioned (immediately after COMPUTE_VISCOSITY + MESH_EXCHANGE(4)). |
| RADIATION before WALL_BC (corrector) | RADIATION reads B1%TMP_G (boundary temperatures) set by WALL_BC. Hard dependency for thermally-thick walls. |

### Revised Pipeline Diagrams

**Predictor (with reordering):**

```
INSERT_PARTICLES → COMPUTE_VISCOSITY → MASS_FINITE_DIFF → DENSITY
  → MESH_EXCHANGE(1)
  → [two-way fork]
     Branch A: VELOCITY_FLUX → PARTICLE_MOMENTUM (266 ops/cell, critical)
     Branch B: WALL_BC (65 ops/cell, hidden behind A)
  → Join → DIV_P1 (full) → DIV_EXCHANGE → DIV_P2 → PRESSURE_SOLVE
  → VEL_PREDICTOR → CHECK_STABILITY → MESH_EXCHANGE(3) → VEL_BC
```

> Note: Three-way fork (VFLUX || WALL_BC || DIV_P1 interior) was investigated and ruled out — DP *= RTRM makes interior/wall split impossible.

**Corrector (with reordering):**

```
COMPUTE_VISCOSITY → MASS_FINITE_DIFF → DENSITY
  → MESH_EXCHANGE(4)
  → [two-way fork]
     Branch A: VELOCITY_FLUX ──────────────────→ PARTICLE_MOMENTUM
     Branch B: COMBUSTION → COND → PME → MOVE ─/
  → EX(7) → WALL_BC → EX(6)
  → [two-way fork]
     Branch C: RADIATION ────────────────→ DIV_P1 QR addition
     Branch D: DIV_P1 non-QR sections ──/
  → DIV_EXCHANGE → DIV_P2
  → [two-way fork]
     Branch E: PRHS computation ──→ FFT solve
     Branch F: Pressure BC setup ─/
  → VEL_CORRECTOR → EX(6) → VEL_BC → DUMP
```

### Impact Summary

| Opportunity | Phase | Est. Benefit | Complexity | Priority |
|-------------|-------|-------------|------------|----------|
| VFLUX \|\| COMBUSTION chain | Corrector | **Moderate** -- 125 ops/cell hidden | Low (no WORK conflict) | 1 |
| DIV_P1 non-QR \|\| RADIATION | Corrector | **High** -- 1299 ops/cell hidden behind RADIATION | Medium (extract QR, scratch dup.) | 2 |
| VFLUX \|\| WALL_BC | Predictor | **Low** -- 65 ops/cell hidden | Low (reuse fork/join pattern) | 3 |
| ~~DIV_P1 interior \|\| WALL_BC~~ | ~~Predictor~~ | ~~331 ops/cell~~ | ~~Medium~~ | ~~RULED OUT~~ |
| Pressure BC \|\| DIV_P2 | Both | **Low** -- BC setup is ~20% of pre-solve | Low (already separate loop) | 4 |

The corrector fork 2 (RADIATION || DIV_P1) is by far the most impactful opportunity, hiding the entire DIV_P1 cost behind RADIATION. The predictor three-way fork was ruled out because `DP *= RTRM` makes an interior/wall split impossible (see Phase 6 investigation).

## Section-Level Cost Analysis and Hedgehog Graph Recommendations

### Reference Scenario

All costs are estimated as **floating-point operations per mesh cell** for a reference configuration:

- Mesh: 64³ = 262,144 cells
- Species: NS = 5 (N_TOTAL_SCALARS)
- Turbulence: CONSMAG (no global filters)
- Combustion: 10% active cells
- Wall cells: ~10% of mesh cells (6 × 64² = 24,576 boundary cells)
- Radiation: NUMBER_RADIATION_ANGLES = 104, ANGLE_INCREMENT = 15
- No CC_IBM, no particles

MFLOPs = ops/cell × 262,144 / 10⁶. Costs labeled with "NS × K" scale linearly with species count.

### Cost Table

#### Pre-Pipeline (Sequential, Already Parallelized via Mesh + K-Block)

Executed once per predictor and once per corrector phase (total: 2 × 1490 = 2980 ops/cell per timestep).

| Section | Ops/cell | MFLOPs | Scaling | Bottleneck |
|---------|----------|--------|---------|------------|
| INSERT_PARTICLES | ~0 | ~0 | per particle | N/A |
| COMPUTE_VISCOSITY: MU_DNS | 15 | 4 | ~const | Mixture viscosity |
| COMPUTE_VISCOSITY: Strain rate | 73 | 19 | ~const | 9 velocity diffs + 6 squares |
| COMPUTE_VISCOSITY: Turbulent MU | 5 | 1 | ~const | CONSMAG: 1 multiply |
| COMPUTE_VISCOSITY: KRES | 12 | 3 | ~const | 3 interpolations + 3 squares |
| **COMPUTE_VISCOSITY total** | **90** | **24** | | |
| MASS_FINITE_DIFF | 400 | 105 | NS × 80 | GET_SCALAR_FACE_VALUE flux limiter |
| DENSITY | 1000 | 262 | NS × 180 | Species flux limiter + density sum |
| **Pre-pipeline total** | **1490** | **391** | | |

#### Predictor Pipeline (After MESH_EXCHANGE(1))

> Three-way fork (VFLUX || WALL_BC || DIV_P1 interior) was investigated and ruled out. `DP *= RTRM` in Phase 5 makes interior/wall split impossible. See [IMPLEMENTATION_PROGRESS.md](IMPLEMENTATION_PROGRESS.md) Phase 6.

| Section | Ops/cell | MFLOPs | Pipeline Branch | Critical Path? |
|---------|----------|--------|-----------------|----------------|
| VFLUX: Vorticity + stress | 44 | 12 | A | **Yes** |
| VFLUX: FVX/FVY/FVZ | 171 | 45 | A | **Yes** |
| VFLUX: Baroclinic corr. | 33 | 9 | A | **Yes** |
| PARTICLE_MOMENTUM | 18 | 5 | A (after VFLUX) | **Yes** |
| **Branch A total** | **266** | **70** | **Critical path** | |
| WALL_BC (all sub-kernels) | 65† | 17 | B | No (hidden behind A) |
| **Branch B total** | **65** | **17** | **Hidden behind A** | |
| DIV_P1 (full) | 1300 | 341 | Sequential (after join) | **Yes** |
| DIV_P2 | 40 | 10 | Sequential | **Yes** |
| PRESSURE_SOLVE | 50 | 13 | Sequential | **Yes** |
| VEL_PREDICTOR | 18 | 5 | Sequential | **Yes** |
| **Predictor pipeline total** | **1739** | **456** | | |

† WALL_BC is 210–1175 ops per wall cell; 65 ops/cell effective assumes ~10% wall fraction with simple surfaces. Thermally-thick walls with SOLID_HEAT_TRANSFER increase this substantially.

**Predictor speedup**: Fork segment 331 → 266 ops/cell = **1.24× fork**. Full pipeline: 1739 → 1674 ops/cell = **1.04×** (WALL_BC is a small fraction of the predictor).

#### Corrector Pipeline (After MESH_EXCHANGE(4))

| Section | Ops/cell | MFLOPs | Pipeline Branch | Critical Path? |
|---------|----------|--------|-----------------|----------------|
| VELOCITY_FLUX | 248 | 65 | A | **Yes** (dominates B) |
| **Branch A total** | **248** | **65** | | |
| COMBUSTION | 115 | 30 | B | No |
| CONDENSATION | 10 | 3 | B | No |
| PART_MASS_ENERGY | varies | - | B | No |
| MOVE_PARTICLES | varies | - | B | No |
| **Branch B total** | **125+** | **33+** | **Hidden behind A** | |
| PARTICLE_MOMENTUM | 18 | 5 | After A+B join | **Yes** |
| WALL_BC | 65 | 17 | Sequential | **Yes** |
| RADIATION | 3500‡ | 917 | C | **Yes** (dominates D) |
| **Branch C total** | **3500** | **917** | | |
| DIV_P1 (excl. QR line) | 1299 | 340 | D | No (hidden behind C) |
| **Branch D total** | **1299** | **340** | **Hidden behind C** | |
| DIV_P1: QR addition | 1 | ~0 | After C+D join | Negligible |
| DIV_P2 | 40 | 10 | Sequential | **Yes** |
| PRESSURE_SOLVE | 50 | 13 | Sequential | **Yes** |
| VEL_CORRECTOR | 18 | 5 | Sequential | **Yes** |
| **Corrector pipeline total** | **5364** | **1406** | | |

‡ RADIATION cost varies enormously: 500 ops/cell (gray gas, few angles) to 15,000 ops/cell (spectral with many bands). 3500 is typical for a multi-band fire scenario.

**Corrector speedup**: Sequential 5364 → Pipelined 3940 ops/cell = **1.36×**

#### Overall Timestep Speedup

| Phase | Sequential | Pipelined | Saved | Speedup |
|-------|-----------|-----------|-------|---------|
| Predictor pipeline | 1739 | 1674 | 65 | 1.04× |
| Corrector pipeline | 5364 | 3940 | 1424 | 1.36× |
| **Combined pipeline** | **7103** | **5614** | **1489** | **1.27×** |
| Full timestep (incl. pre-pipeline 2×1490) | 10083 | 8594 | 1489 | **1.17×** |

### Critical Path Analysis

The critical path determines execution time. Sections not on the critical path run "for free" behind the bottleneck.

**Predictor critical path** (1674 ops/cell):

```
VFLUX(248) + PART_MOM(18) → DIV_P1(1300) → zone_sums(15)
  → DIV_EXCHANGE → DIV_P2(40) → PRES(50) → VEL_PRED(18)
```

Hidden (concurrent with VFLUX+PART_MOM):
- WALL_BC (65 ops) — completes at 24% of Branch A time

> Note: Three-way fork (hiding DIV_P1 interior behind WALL_BC) was ruled out — `DP *= RTRM` prevents interior/wall split.

**Corrector critical path** (3940 ops/cell):

```
VFLUX(248) → PART_MOM(18) → EX(7) → WALL_BC(65) → EX(6)
  → RADIATION(3500) → EX(2) → QR_add(1)
  → DIV_EXCHANGE → DIV_P2(40) → PRES(50) → VEL_CORR(18)
```

Hidden:
- COMBUSTION+CONDENSATION (125 ops) — behind VFLUX at Fork 1
- DIV_P1 main (1299 ops) — behind RADIATION at Fork 2 (37% of RADIATION time)

### Sensitivity to RADIATION Cost

The corrector pipelining benefit depends on the relative cost of RADIATION vs DIV_P1:

| RADIATION (ops/cell) | Scenario | Sequential | Pipelined | Speedup |
|---------------------|----------|-----------|-----------|---------|
| 500 | Gray gas | 2494 | 1869 | 1.33× |
| 2000 | Moderate | 3994 | 2869 | 1.39× |
| 3500 | Typical fire | 5364 | 3940 | 1.36× |
| 5000 | Spectral | 6864 | 5540 | 1.24× |
| 15000 | Many bands | 16864 | 15540 | 1.09× |

When RADIATION cost is comparable to DIV_P1 (1000–3000 ops/cell), pipelining provides the best relative speedup: whichever routine is shorter runs entirely hidden behind the longer one. When RADIATION dominates (>5000), DIV_P1 is always hidden but the relative savings shrink since RADIATION is the bottleneck regardless.

### Hedgehog Graph Recommendations

Ranked by impact-to-complexity ratio. See [fds_section_pipeline.dot](fds_section_pipeline.dot) for the visual structure.

#### Tier 1 — High Impact, Implement First

**1.1 Corrector: VFLUX || COMBUSTION chain (after EX4)**

Fork into two branches after MESH_EXCHANGE(4):
- Branch A: VelocityFluxTask (248 ops/cell)
- Branch B: CombustionTask → CondensationTask → PartMassEnergyTask → MoveParticlesTask (125+ ops/cell)

Join at: ParticleMomentumTask (needs FVX from A, FVX_D from B).

Implementation: Type-based fork. Emit `VelocityFluxWork` and `CombustionChainWork` from a fork state. Collector waits for both result types per mesh.

WORK conflict: None. COMBUSTION and CONDENSATION use no WORK arrays. PME starts after VFLUX completes (short branch). With scratch duplication, even timing-dependent overlaps are safe.

Savings: 125 ops/cell minimum. With particles, Branch B grows longer, making the overlap increasingly "free."

**1.2 Corrector: RADIATION || DIV_P1 (after EX6)**

Fork into two branches after MESH_EXCHANGE(6):
- Branch C: RadiationTask (500–15000 ops/cell)
- Branch D: DivP1NonQRTask (1299 ops/cell) — everything except the single QR addition line

Join at: DivP1QRAdditionTask (applies `DP += QR`, 1 op/cell).

Implementation: Split DIVERGENCE_PART_1_KERNEL into two Fortran subroutines:
- `DIV_P1_MAIN_KERNEL`: All phases (2A, 2B, 3, 4, 5A–G, 6) with QR set to zero
- `DIV_P1_QR_ADDITION`: Single loop: `DP(I,J,K) = DP(I,J,K) + QR(I,J,K)`

WORK conflict: RADIATION uses WORK1-9, DIV_P1 uses WORK1-7,9. **Requires per-branch scratch duplication** (second set of WORK arrays on each mesh, ~4.5 MB per 64³ mesh).

Savings: min(RADIATION_cost, 1299) ops/cell. For typical fire: 1299 ops/cell saved — the entire DIV_P1 runs hidden behind RADIATION.

#### Tier 2 — Moderate Impact, Implement Next

**2.1 Predictor: Two-way fork (VFLUX || WALL_BC)**

Fork into two branches after MESH_EXCHANGE(1):
- Branch A: VelocityFluxTask → ParticleMomentumTask (266 ops/cell, critical path)
- Branch B: WallBCTask (65 ops/cell, hidden behind A)

Join at: DIV_P1 (full computation, needs B1 props from WALL_BC).

No WORK conflict (WALL_BC uses no WORK arrays). No Fortran kernel changes needed. Reuses fork/join pattern from Tier 1.1.

Savings: 65 ops/cell (WALL_BC hidden behind VFLUX+PART_MOM).

> **Three-way fork (VFLUX || WALL_BC || DIV_P1_interior) was ruled out.** DIV_P1's interior I,J,K loops include boundary face positions that are read into DP. The multiplicative step `DP *= RTRM` in Phase 5 makes any boundary-face error irreversible. See [IMPLEMENTATION_PROGRESS.md](IMPLEMENTATION_PROGRESS.md) Phase 6.

#### Tier 3 — Low Impact, Deferred

**3.1 Pressure BC || DIV_P2** (both phases): Saves ~3 ops/cell. Not worth the graph complexity.

**3.2 Intra-DIV_P1 Groups A||B||C**: Phases 2A+2B (200), 3+4 (95), 5A (105) can run in parallel. Saves ~200 ops/cell within DIV_P1's first half. High complexity (3 separate WORK pools, DP accumulation synchronization). Diminishing returns since DIV_P1 is already partially hidden.

**3.3 VFLUX FVX||FVY||FVZ**: Saves ~114 ops/cell within VFLUX. But VFLUX is already fully hidden in the predictor (behind DIV_P1). In the corrector, VFLUX is on the critical path but not the bottleneck.

#### Out of Scope for Hedgehog (Use OpenMP Instead)

**RADIATION angle parallelism**: Each of the 100–500 angles sweeps independently. Potential 10–50× within RADIATION. Requires per-thread private arrays (IL, IL_UP). This is loop-level parallelism — better suited to OpenMP DO PARALLEL than Hedgehog graph structure.

**COMBUSTION cell-level parallelism**: Embarrassingly parallel per active cell. Fine-grained — better suited to OpenMP or GPU.

### Implementation Roadmap

```
Phase 1: Corrector VFLUX || COMBUSTION fork (Tier 1.1)
  - No Fortran kernel changes needed
  - Add fork state after EX(4), join state before PART_MOM
  - Estimated effort: low
  - Expected speedup: ~2% of corrector, higher with particles

Phase 2: Corrector RADIATION || DIV_P1 fork (Tier 1.2)
  - Fortran: split DIV_P1 into main + QR-addition kernels
  - Add per-branch WORK array pools to MESH_TYPE (one-time, 4.5 MB/mesh)
  - Add fork state after EX(6), join state before QR addition
  - Estimated effort: medium
  - Expected speedup: ~30% of corrector pipeline

Phase 3: Predictor two-way fork (Tier 2.1)
  - No Fortran kernel changes needed
  - Reuse fork/join pattern from Phase 1
  - Add fork state after EX(1), join before DIV_P1
  - Estimated effort: low
  - Expected speedup: ~4% of predictor pipeline (65 ops/cell saved)

Phase 4 (optional): Intra-routine optimizations (Tier 3)
  - Only pursue after Phases 1-3 are validated
  - Profile to identify actual bottlenecks vs estimates
```

### DIV_P1 Phase Breakdown

Detailed cost breakdown of DIVERGENCE_PART_1_KERNEL phases, the most expensive routine:

| Phase | Interior Ops/cell | Wall Corr. Ops/cell | Description | WORK Usage |
|-------|-------------------|---------------------|-------------|------------|
| 1 | 0 | 0 | Setup (zero DP, aliases) | — |
| 2A | 100 | 30 | Species diffusion fluxes (RHO_D_DZDX/Y/Z) | WORK1-3 |
| 2B | 55 | 15 | Diffusive heat flux (H_RHO_D_DZDX/Y/Z → DP) | WORK4-6 |
| 3 | 35 | 0 | Specific heat (CP, R_H_G) | WORK7 |
| 4 | 50 | 10 | Thermal conductivity + divergence (→ DP + **QR**) | WORK1-3 (reused) |
| 5A | 80 | 25 | Enthalpy advection (→ DP) | WORK1-3 (reused) |
| 5B | 5 | 0 | RTRM = 1/(rho×CP×TMP); DP *= RTRM | — |
| 5C | 670 | 0 | Species face flux limiter (GET_SCALAR_FACE_VALUE) | SWORK1-3 |
| 5D | 60 | 40 | Species divergence (→ DP) | WORK9 |
| 5E-G | 10 | 0 | Source terms (reactions, stratification, MMS) | — |
| 6 | 10 | 5 | Pressure zone sums (DSUM, PSUM, USUM) | — |
| **Total** | **~1075** | **~125** | | |

Phase 5C (species advection flux limiter) is the single most expensive section at ~670 ops/cell for NS=5. It dominates DIV_P1 cost at 50%+. GET_SCALAR_FACE_VALUE computes 40–60 FLOPs per face value for flux-limited interpolation.

The interior/wall correction split is clean: bulk I,J,K loops compute on interior cells using mesh arrays available before the pipeline fork. Wall correction loops iterate over wall cells and read B1 boundary properties (TMP_F, RHO_F, ZZ_F, U_NORMAL_S, RHO_D_F) set by WALL_BC.

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
