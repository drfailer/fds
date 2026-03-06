# FDS Architecture Analysis

This document provides a comprehensive analysis of the Fire Dynamics Simulator (FDS)
codebase architecture, including module dependencies, execution flow, data structures,
thread safety assessment, and inter-process communication patterns.

---

## Table of Contents

1. [Module Dependency Graph](#1-module-dependency-graph)
2. [Execution Flow](#2-execution-flow)
3. [Data Types and Data Flow](#3-data-types-and-data-flow)
4. [Intra-Node Parallelism (Thread Safety)](#4-intra-node-parallelism-thread-safety)
5. [Inter-Node Parallelism (MPI Communication)](#5-inter-node-parallelism-mpi-communication)

---

## 1. Module Dependency Graph

### 1.1 Module Inventory

FDS is organized into ~50 Fortran modules across the following categories:

| Category | Modules | Files |
|----------|---------|-------|
| **Core Infrastructure** | PRECISION_PARAMETERS, GLOBAL_CONSTANTS, TYPES, MESH_VARIABLES, OUTPUT_DATA | prec.f90, cons.f90, type.f90, mesh.f90, data.f90 |
| **Numerical Solvers** | PRES, PRES_KERNELS, POIS, VELO, VELO_KERNELS, DIVG, DIVG_KERNELS, MASS, MASS_KERNELS | pres.f90, pres_kernels.f90, pois.f90, velo.f90, velo_kernels.f90, divg.f90, divg_kernels.f90, mass.f90, mass_kernels.f90 |
| **Physics** | FIRE, FIRE_KERNELS, RAD, TURBULENCE, TURB_KERNELS, WALL_ROUTINES, WALL_KERNELS, PART, CVODE_INTERFACE, HVAC_ROUTINES, SOOT_ROUTINES, VEGE | fire.f90, fire_kernels.f90, radi.f90, turb.f90, turb_kernels.f90, wall.f90, wall_kernels.f90, part.f90, chem.f90, hvac.f90, soot.f90, vege.f90 |
| **I/O & Setup** | READ_INPUT, INIT, DUMP, COMP_FUNCTIONS, DEVICE_VARIABLES, CONTROL_FUNCTIONS | read.f90, init.f90, dump.f90, func.f90, devc.f90, ctrl.f90 |
| **Complex Geometry** | CC_SCALARS_DATA, CC_SCALARS, CC_INIT, CC_DIVERGENCE, CC_DIVERGENCE_KERNELS, CC_VELOCITY, CC_VELOCITY_KERNELS, CC_PRESSURE, CC_PRESSURE_KERNELS, CC_DENSITY_MOD, CC_EXCHANGE, CC_VERIFICATION, COMPLEX_GEOMETRY | ccib/*.f90, geom.f90 |
| **Utilities** | MATH_FUNCTIONS, PHYSICAL_FUNCTIONS, MEMORY_FUNCTIONS, MISC_FUNCTIONS, MANUFACTURED_SOLUTIONS, TRAN, THERMO_PROPS, RADCAL_VAR, BOXTETRA_ROUTINES, ISOSMOKE | Various |
| **Hedgehog** | FDS_DRIVER, FDS_C_INTERFACE | hedgehog/fds_driver.f90, hedgehog/fds_c_interface.f90 |

### 1.2 Orchestration / Kernel Pattern

Ten modules follow a split pattern where the **orchestration** module manages mesh
pointers and control flow, while the **kernel** module contains pure computation
routines taking `TYPE(MESH_TYPE), INTENT(INOUT)` as explicit argument:

| Orchestration | Kernel | Domain |
|---------------|--------|--------|
| VELO | VELO_KERNELS | Velocity prediction/correction |
| DIVG | DIVG_KERNELS | Divergence calculation |
| MASS | MASS_KERNELS | Density/species transport |
| PRES | PRES_KERNELS | Pressure solver RHS & FFT |
| FIRE | FIRE_KERNELS | Combustion |
| TURBULENCE | TURB_KERNELS | Turbulence modeling |
| WALL_ROUTINES | WALL_KERNELS | Wall boundary conditions |
| CC_DIVERGENCE | CC_DIVERGENCE_KERNELS | Cut-cell divergence |
| CC_VELOCITY | CC_VELOCITY_KERNELS | Cut-cell velocity |
| CC_PRESSURE | CC_PRESSURE_KERNELS | Cut-cell pressure |

### 1.3 Dependency Graph (DOT)

See `graphs/module_dependencies.dot` for the full graph. Render with:
```bash
dot -Tsvg graphs/module_dependencies.dot -o graphs/module_dependencies.svg
```

The dependency layers are:

```
Layer 0 (Foundation):   PRECISION_PARAMETERS
                              |
Layer 1 (Constants):    GLOBAL_CONSTANTS, MPI_F08
                              |
Layer 2 (Types):        TYPES, MKL_PARDISO, HYPRE_INTERFACE
                              |
Layer 3 (Mesh):         MESH_VARIABLES (defines MESH_TYPE)
                              |
Layer 4 (Utilities):    MATH_FUNCTIONS, PHYSICAL_FUNCTIONS, COMP_FUNCTIONS,
                        OUTPUT_DATA, MEMORY_FUNCTIONS, TRAN
                              |
Layer 5 (Kernels):      VELO_KERNELS, DIVG_KERNELS, MASS_KERNELS,
                        PRES_KERNELS, FIRE_KERNELS, TURB_KERNELS,
                        WALL_KERNELS, CC_*_KERNELS
                              |
Layer 6 (Orchestration): VELO, DIVG, MASS, PRES, FIRE, TURBULENCE,
                         WALL_ROUTINES, CC_DIVERGENCE, CC_VELOCITY,
                         CC_PRESSURE
                              |
Layer 7 (Physics):      RAD, PART, HVAC_ROUTINES, SOOT_ROUTINES,
                        VEGE, CC_SCALARS, CC_DENSITY_MOD
                              |
Layer 8 (I/O):          READ_INPUT, INIT, DUMP
                              |
Layer 9 (Driver):       main.f90 / FDS_DRIVER
```

Key dependency rules:
- **Kernel modules** depend only on: PRECISION_PARAMETERS, TYPES, GLOBAL_CONSTANTS, MESH_VARIABLES, utility modules
- **Orchestration modules** depend on their kernel + MESH_POINTERS + COMP_FUNCTIONS
- **Physics modules** may depend on other physics modules (e.g., DIVG_KERNELS uses TURB_KERNELS)
- **I/O modules** depend on everything (broadest dependency set)

---

## 2. Execution Flow

### 2.1 High-Level Program Structure

```
PROGRAM FDS
  |
  +-- [Initialization Phase]     ~500 lines
  |     MPI setup, input parsing, mesh allocation, atmosphere,
  |     wall arrays, complex geometry, radiation, particles
  |
  +-- [Main Time-Stepping Loop]  ~500 lines
  |     |
  |     +-- PREDICTOR step
  |     |     Insert particles, compute viscosity, mass transport,
  |     |     velocity flux, wall BC, divergence, pressure solve,
  |     |     velocity prediction, CFL check
  |     |
  |     +-- CORRECTOR step
  |     |     Mass transport correction, combustion, particle transfer,
  |     |     wall BC, radiation, divergence, pressure solve,
  |     |     velocity correction, output/diagnostics
  |     |
  |     +-- [Output and diagnostics]
  |
  +-- [Finalization]             ~20 lines
        Solver cleanup, MPI finalize
```

### 2.2 Execution Flow Graph (DOT)

See `graphs/execution_flow.dot` for the complete graph. Key aspects:

**Within each mesh loop**, the computation order is:

```
PREDICTOR:
  INSERT_ALL_PARTICLES --> COMPUTE_VISCOSITY --> MASS_FINITE_DIFFERENCES_NEW
  --> DENSITY --> [MESH_EXCHANGE(1)] --> VELOCITY_FLUX --> WALL_BC
  --> DIVERGENCE_PART_1 --> [EXCHANGE_DIVERGENCE_INFO]
  --> DIVERGENCE_PART_2 --> PRESSURE_ITERATION_SCHEME
  --> VELOCITY_PREDICTOR --> CHECK_STABILITY
  --> [MESH_EXCHANGE(3)] --> MATCH_VELOCITY --> VELOCITY_BC

CORRECTOR:
  COMPUTE_VISCOSITY --> MASS_FINITE_DIFFERENCES_NEW --> DENSITY
  --> [MESH_EXCHANGE(4)] --> VELOCITY_FLUX
  --> COMBUSTION_LOAD_BALANCED (global, not per-mesh)
  --> PARTICLE_MASS_ENERGY_TRANSFER --> MOVE_PARTICLES
  --> [MESH_EXCHANGE(11)] --> WALL_BC --> COMPUTE_RADIATION
  --> [MESH_EXCHANGE(2)] --> DIVERGENCE_PART_1
  --> [EXCHANGE_DIVERGENCE_INFO] --> DIVERGENCE_PART_2
  --> PRESSURE_ITERATION_SCHEME --> VELOCITY_CORRECTOR
  --> [MESH_EXCHANGE(6)] --> MATCH_VELOCITY --> VELOCITY_BC
  --> UPDATE_GLOBAL_OUTPUTS --> DUMP_MESH_OUTPUTS
```

### 2.3 Pressure Iteration Detail

The pressure solver is the most complex inner loop:

```
PRESSURE_ITERATION_LOOP:
  |
  +-- BAROCLINIC_CORRECTION(NM)     [per-mesh]
  +-- CC_NO_FLUX(NM)                [per-mesh, if CC_IBM]
  +-- MESH_EXCHANGE(5)              [MPI sync]
  +-- MATCH_VELOCITY_FLUX(NM)       [per-mesh]
  +-- NO_FLUX(NM)                   [per-mesh]
  +-- PRESSURE_SOLVER_COMPUTE_RHS   [per-mesh]
  +-- SELECT CASE(PRES_FLAG):
  |     FFT_FLAG  --> PRESSURE_SOLVER_FFT(NM)     [per-mesh]
  |     GLMAT_FLAG --> GLMAT_SOLVER(T,DT)          [global, all meshes]
  |     ULMAT_FLAG --> ULMAT_SOLVER(NM)            [per-mesh]
  +-- PRESSURE_SOLVER_CHECK_RESIDUALS(NM)  [per-mesh]
  +-- if ITERATE_PRESSURE:
        MESH_EXCHANGE(5)
        COMPUTE_VELOCITY_ERROR(NM)          [per-mesh]
        MPI_ALLGATHERV (velocity errors)    [global sync]
        Check convergence --> EXIT or CYCLE
```

### 2.4 Time Step Adaptation

The predictor wraps the density-through-velocity computation in a
`CHANGE_TIME_STEP_LOOP`. If CFL check fails after `VELOCITY_PREDICTOR`:
1. `CHANGE_TIME_STEP_INDEX = -1` for the unstable mesh
2. `MPI_ALLGATHERV` synchronizes the flag across processes
3. DT is reduced, and the entire predictor restarts from DENSITY

---

## 3. Data Types and Data Flow

### 3.1 MESH_TYPE - The Central Data Structure

`MESH_TYPE` (defined in mesh.f90) is the primary per-mesh container. Each MPI
process holds a `MESHES(1:NMESHES)` array, but only processes meshes
`LOWER_MESH_INDEX` through `UPPER_MESH_INDEX`.

#### Key 3D Field Arrays (I,J,K indexed)

| Array | Physical Quantity | Units | Updated By |
|-------|------------------|-------|------------|
| `U, V, W` | Velocity components (face-centered) | m/s | VELO (predictor/corrector) |
| `US, VS, WS` | Estimated velocity at next time step | m/s | VELO (predictor) |
| `H, HS` | Stagnation pressure: p'/rho + \|u\|^2/2 | m^2/s^2 | PRES |
| `RHO, RHOS` | Gas density (current/estimated) | kg/m^3 | MASS |
| `TMP` | Gas temperature | K | MASS (from EOS) |
| `D, DS, DDDT` | Divergence, estimated, time derivative | 1/s | DIVG |
| `MU` | Dynamic (turbulent) viscosity | Pa.s | TURB |
| `Q` | Heat release rate per unit volume | W/m^3 | FIRE |
| `QR, QR_W` | Radiation source term | W/m^3 | RAD |
| `KRES` | Resolved kinetic energy | m^2/s^2 | VELO |
| `PRHS` | Poisson equation RHS | 1/s^2 | PRES |
| `FVX, FVY, FVZ` | Momentum flux terms | m/s^2 | VELO |

#### Key 4D Arrays (I,J,K,N indexed)

| Array | Physical Quantity | Updated By |
|-------|------------------|------------|
| `ZZ(I,J,K,N), ZZS` | Species mass fractions (current/estimated) | MASS |
| `FX, FY, FZ(I,J,K,N)` | Convective species fluxes | MASS |
| `DEL_RHO_D_DEL_Z(I,J,K,N)` | Diffusive species flux divergence | MASS |
| `REAC_SOURCE_TERM(I,J,K,N)` | Chemical reaction source terms | FIRE |
| `M_DOT_PPP(I,J,K,N)` | Particle mass source per unit volume | PART |

#### Boundary Data (Indexed by wall cell IW)

| Component Array | Purpose |
|----------------|---------|
| `WALL(IW)` | Wall cell metadata (indices into component arrays) |
| `BOUNDARY_COORD(BC)` | Ghost/gas cell indices, normal vector, orientation |
| `BOUNDARY_ONE_D(OD)` | 1-D solid conduction/pyrolysis solver state |
| `BOUNDARY_PROP1(B1)` | Surface temperatures, heat/mass transfer coefficients |
| `BOUNDARY_PROP2(B2)` | Droplet accumulation, wall model (u_tau, y+) |
| `BOUNDARY_RADIA(BR)` | Angular radiation intensities at boundaries |

#### Neighbor Data

| Component | Purpose |
|-----------|---------|
| `OMESH(NOM)` | Ghost cell data from neighboring mesh NOM |
| `EXTERNAL_WALL(IW)` | Inter-mesh boundary cell data |

### 3.2 Global (Not Per-Mesh) Data Structures

| Type | Storage | Purpose |
|------|---------|---------|
| `SPECIES_MIXTURE_TYPE` | `SPECIES_MIXTURE(1:N_SMIX)` | Lumped species definitions |
| `REACTION_TYPE` | `REACTION(1:N_REACTIONS)` | Chemical reaction parameters |
| `SURFACE_TYPE` | `SURFACE(0:N_SURF)` | Boundary condition specifications |
| `LAGRANGIAN_PARTICLE_CLASS_TYPE` | `LAGRANGIAN_PARTICLE_CLASS(1:N_PART)` | Particle class definitions |
| `RAMPS_TYPE` | `RAMPS(1:N_RAMP)` | Time/temperature-dependent ramp functions |
| `P_ZONE_TYPE` | `P_ZONE(1:N_ZONE)` | Pressure zone definitions |

These are **read-only during time stepping** (set during initialization).

### 3.3 Data Flow Through One Time Step

See `graphs/data_flow.dot` for the complete graph. The primary data flow is:

```
                    PREDICTOR                              CORRECTOR
                    =========                              =========

ZZ,RHO -----> MASS_FINITE_DIFFERENCES ------+
  |           (computes FX,FY,FZ,           |
  |            DEL_RHO_D_DEL_Z)             |
  |                                         |
  +-------> DENSITY ---------> ZZS, RHOS    +-------> DENSITY --------> ZZ, RHO
              (transport eq)     |                      (correction)      |
                                 |                                        |
U,V,W -----> VELOCITY_FLUX ---> FVX,FVY,FVZ                             |
              (advection +       |                                        |
               diffusion)        |                                        |
                                 v                                        |
ZZS,RHOS --> WALL_BC ---------> Surface fluxes                           |
              (T, species BC)    |                                        |
                                 |                                        |
              DIVERGENCE_PART_1 -> D, DS                                  |
              (heat release,      |                                       |
               mass sources)      |                                       |
                                  v                                       |
              DIVERGENCE_PART_2 -> DDDT                                   |
                                  |                                       |
              PRESSURE_SOLVER --> H, HS                                   |
              (Poisson eq)        |                                       |
                                  v                                       |
              VELOCITY_PREDICTOR -> US, VS, WS                            |
              (U* = U - dt*dH/dx)  |                                      |
                                   |                                      |
              CHECK_STABILITY ---> DT adjustment                          |
                                                                          |
                                                    COMBUSTION ----------> Q, REAC_SOURCE_TERM
                                                    (chemistry ODE)        |
                                                                           |
                                                    PARTICLE_TRANSFER ---> M_DOT_PPP, FVX_D
                                                    (droplet evaporation)  |
                                                                           |
                                                    RADIATION ------------> QR, QR_W
                                                    (RTE solver)           |
                                                                           v
                                                    DIVERGENCE ----------> D, DDDT
                                                                           |
                                                    PRESSURE_SOLVER -----> H
                                                                           |
                                                    VELOCITY_CORRECTOR --> U, V, W (final)
```

### 3.4 Key Data Transformation Chain

The fundamental equation solved is the low-Mach number Navier-Stokes:

1. **Species transport** (MASS): `d(rho*Z)/dt + div(rho*Z*u) = div(rho*D*grad(Z)) + sources`
   - Input: `ZZ, RHO, U, V, W` --> Output: `ZZS, RHOS` (predicted) or `ZZ, RHO` (corrected)

2. **Equation of state** (embedded in MASS): `rho = P_0 * W_mix / (R * T)`
   - Input: `ZZS, RHOS, PBAR` --> Output: `TMP` (temperature)

3. **Divergence** (DIVG): `div(u) = (1/rho) * [sum of heat/mass sources]`
   - Input: `Q, QR, M_DOT_PPP, ZZ, TMP, PBAR` --> Output: `D, DDDT`

4. **Poisson equation** (PRES): `Laplacian(H) = -dD/dt - div(F)`
   - Input: `D, DDDT, FVX, FVY, FVZ` --> Output: `H` (pressure head)

5. **Velocity update** (VELO): `u^{n+1} = u* - dt * grad(H) / rho`
   - Input: `FVX, FVY, FVZ, H, RHO` --> Output: `U, V, W`

---

## 4. Intra-Node Parallelism (Thread Safety)

The goal is to process multiple meshes concurrently within a single MPI process
(task-based parallelism, NOT OpenMP). This section identifies all barriers to
concurrent per-mesh execution.

### 4.1 Threat Categories

#### CRITICAL: POINT_TO_MESH Pattern

`POINT_TO_MESH(NM)` (mesh.f90:356-486) sets ~200 module-level pointer aliases
to components of `MESHES(NM)`. This is **fundamentally thread-unsafe**: if two
threads call `POINT_TO_MESH` for different meshes, they corrupt each other's
pointers.

**Status**: 14 modules still use POINT_TO_MESH:

| Module | File | Uses |
|--------|------|------|
| PRES | pres.f90 | Orchestration calls |
| VELO | velo.f90 | Orchestration calls |
| DIVG | divg.f90 | Orchestration calls |
| MASS | mass.f90 | Orchestration calls |
| FIRE | fire.f90 | Combustion orchestration |
| WALL_ROUTINES | wall.f90 | Wall BC orchestration |
| TURBULENCE | turb.f90 | Turbulence verification |
| PART | part.f90 | Particle tracking |
| RAD | radi.f90 | Radiation transport |
| HVAC_ROUTINES | hvac.f90 | HVAC solver |
| SOOT_ROUTINES | soot.f90 | Soot agglomeration |
| VEGE | vege.f90 | Level-set fire spread |
| DUMP | dump.f90 | Output routines |
| READ_INPUT | read.f90 | Input parsing (init only) |

**Already safe** (kernel modules with explicit `M%` access):
VELO_KERNELS, DIVG_KERNELS, MASS_KERNELS, PRES_KERNELS, FIRE_KERNELS,
TURB_KERNELS, WALL_KERNELS, CC_DIVERGENCE_KERNELS, CC_VELOCITY_KERNELS,
CC_PRESSURE_KERNELS

**Mitigation**: Continue kernel extraction. The orchestration modules (VELO,
DIVG, MASS, PRES, FIRE, WALL_ROUTINES) are thin wrappers that call kernels,
so refactoring them is tractable. PART and RAD are the largest remaining
challenges.

#### CRITICAL: Global Accumulators (OUTPUT_DATA)

```fortran
! data.f90
REAL(EB), ALLOCATABLE :: Q_DOT(:)      ! Heat release rate accumulator
REAL(EB), ALLOCATABLE :: M_DOT(:)      ! Mass loss rate accumulator
REAL(EB), ALLOCATABLE :: Q_DOT_SUM(:)  ! Time-integrated HRR
REAL(EB), ALLOCATABLE :: M_DOT_SUM(:)  ! Time-integrated mass
```

These are **zeroed at the start of the corrector** (`Q_DOT=0; M_DOT=0` in
main.f90) and then **accumulated across all meshes** during DIVERGENCE_PART_1
and other routines. Concurrent accumulation = data race.

**Mitigation**: Per-mesh accumulation arrays, merged after the mesh loop.

#### CRITICAL: SAVE Variables in Computation Routines

| Module | Variables | Risk |
|--------|-----------|------|
| pres.f90 | `CYL_FCT`, `ILO_CELL..KHI_FACE` (loop bounds) | Written per-mesh |
| geom.f90 | `CC_NEDGECROSS`, `CC_NCUTCELL`, `T_CC_USED`, `ILO_CELL..KHI_FACE` | Written per-mesh, accumulated |
| ccib_data.f90 | `ILO_CELL..KHI_FACE`, `NXB, NYB, NZB` | Written per-mesh |
| pois.f90 | `POIS_WORK` type (FFT state) | Modified during solve |
| fire.f90 | `T_CHEM_ODE`, `COMBUSTION_INIT` | Modified during combustion |
| soot.f90 | `BIN_S, BIN_M`, other arrays | Modified per-mesh |

**Mitigation**:
- Loop bounds: compute locally or pass as arguments
- `CYL_FCT`: pass as argument
- `POIS_WORK`: already addressed (thread-local POIS_WORK type)
- Accumulators: per-mesh with merge

#### MODERATE: Timing Accumulator

```fortran
! cons.f90
REAL(EB), ALLOCATABLE :: T_USED(:)  ! e.g., T_USED(5) = T_USED(5) + elapsed
```

Incremented by every major routine. Race condition under concurrency.

**Mitigation**: Thread-local timing arrays, merged at barrier points.

### 4.2 Thread Safety Summary by Routine

For the per-mesh computation routines called in the main loop:

| Routine | Module | Kernel Safe? | Orchestration Safe? | Blockers |
|---------|--------|:---:|:---:|-----------|
| DENSITY | MASS | Yes | No | POINT_TO_MESH |
| MASS_FINITE_DIFFERENCES_NEW | MASS | Yes | No | POINT_TO_MESH |
| COMPUTE_VISCOSITY | TURB | Yes | No | POINT_TO_MESH |
| VELOCITY_FLUX | VELO | Yes | No | POINT_TO_MESH |
| WALL_BC | WALL | Partial | No | POINT_TO_MESH, module pointers |
| DIVERGENCE_PART_1 | DIVG | Yes | No | POINT_TO_MESH, Q_DOT/M_DOT |
| DIVERGENCE_PART_2 | DIVG | Yes | No | POINT_TO_MESH |
| PRESSURE_SOLVER_FFT | PRES | Yes | No | POINT_TO_MESH, CYL_FCT |
| VELOCITY_PREDICTOR | VELO | Yes | No | POINT_TO_MESH |
| VELOCITY_CORRECTOR | VELO | Yes | No | POINT_TO_MESH |
| INSERT_ALL_PARTICLES | PART | N/A | No | POINT_TO_MESH |
| MOVE_PARTICLES | PART | N/A | No | POINT_TO_MESH |
| COMPUTE_RADIATION | RAD | N/A | No | POINT_TO_MESH, SAVE arrays |
| COMBUSTION_LOAD_BALANCED | FIRE | Partial | No | Global load balancing |
| HVAC_CALC | HVAC | N/A | N/A | Runs on rank 0 only |

### 4.3 Parallelization Roadmap

**Phase 1** (Low-hanging fruit): Refactor the 6 orchestration wrappers (VELO,
DIVG, MASS, PRES, FIRE, WALL_ROUTINES) to pass `MESHES(NM)` directly to
kernels without calling POINT_TO_MESH. These are thin wrappers.

**Phase 2** (Medium effort): Address global accumulators (Q_DOT, M_DOT, T_USED)
with per-mesh scratch arrays and a merge phase.

**Phase 3** (High effort): Refactor PART and RAD to kernel pattern. These are
large modules with deep POINT_TO_MESH usage.

**Phase 4** (Integration): Add task-based scheduling (Hedgehog) to run
per-mesh computations concurrently, with synchronization at MESH_EXCHANGE
boundaries.

---

## 5. Inter-Node Parallelism (MPI Communication)

### 5.1 Process-Mesh Mapping

Each MPI process owns a contiguous range of mesh indices:

```
Process 0: MESHES(1) .. MESHES(K1)
Process 1: MESHES(K1+1) .. MESHES(K2)
...
Process P-1: MESHES(KP-1+1) .. MESHES(NMESHES)
```

Key variables:
- `MY_RANK` - Current process ID (0 to N_MPI_PROCESSES-1)
- `PROCESS(NM)` - Which process owns mesh NM
- `LOWER_MESH_INDEX`, `UPPER_MESH_INDEX` - Mesh range for current process

### 5.2 MESH_EXCHANGE Codes

All inter-mesh data transfer goes through `MESH_EXCHANGE(CODE)` in main.f90.
The CODE parameter specifies what data package to exchange:

| CODE | Data Exchanged | When | Type |
|------|---------------|------|------|
| 0 | Setup persistent MPI requests | Init only | Setup |
| 1 | RHO, MU, KRES, D, Q, ZZ (species/density) | Predictor | Persistent |
| 2 | Radiation intensity (IL_S, IL_R) | After radiation solve | Persistent |
| 3 | H, U, V, W (pressure/velocity, predictor) | End of predictor | Persistent |
| 4 | RHO, MU, D, ZZ (species/density, corrector) | Corrector | Persistent |
| 5 | FVX, FVY, FVZ, H (momentum flux + pressure) | Pressure iteration | Persistent |
| 6 | H, U, V, W + BACK_WALL data | End of corrector | Persistent + blocking |
| 7 | Orphaned particle counts | Before particle exchange | Persistent |
| 8 | Wall cell counts for neighbor setup | Wall init | Blocking |
| 9 | Wall cell indices and SURF_INDEX | Wall init | Blocking |
| 10 | Wall cell data (reals, ints, logicals) | Wall init | Persistent |
| 11 | Particle arrays (REALS, INTEGERS, LOGICALS) | After particle movement | Blocking |
| 14 | Level-set values (PHI_LS, U_LS, V_LS, Z_LS) | If LEVEL_SET_MODE>0 | Persistent |
| 15-18 | Obstruction mass exchange sequence | If EXCHANGE_OBST_MASS | Mixed |
| 19 | Buffer dimension sizes for walls | Wall init | Blocking |

### 5.3 Communication Patterns

**Persistent (non-blocking, reusable)**: Codes 1-6, 10, 14, 15. Initialized
once with `MPI_SEND_INIT`/`MPI_RECV_INIT`, then started with `MPI_STARTALL`
each time step. This is optimal for repeated exchanges.

**Non-persistent (one-time or irregular)**: Codes 8, 9, 11, 16, 18, 19. Use
`MPI_ISEND`/`MPI_IRECV` for irregular data patterns (particles changing mesh).

**Intra-process optimization**: When sender and receiver are on the same MPI
process (`RNODE == SNODE`), data is copied directly without MPI calls.

### 5.4 Global Reductions

Beyond mesh-to-mesh exchange, FDS uses `MPI_ALLREDUCE` for global synchronization:

| Quantity | Operation | When | File |
|----------|-----------|------|------|
| `CHANGE_TIME_STEP_INDEX` | ALLGATHERV | After CFL check | main.f90 |
| `VELOCITY_ERROR_MAX` | ALLGATHERV | Pressure iteration | main.f90 |
| `PRESSURE_ERROR_MAX` | ALLGATHERV | Pressure iteration | main.f90 |
| `STOP_STATUS` | MAX | End of each step | main.f90 |
| `DSUM, PSUM, USUM` | SUM | Zone pressure balance | main.f90 |
| `CONNECTED_ZONES` | MAX | Zone connectivity | main.f90 |
| `EXCHANGE_INSERTED_PARTICLES` | LOR | After particle insertion | main.f90 |
| `EXCHANGE_OBST_MASS` | LOR | Obstruction events | main.f90 |
| `RAD_Q_SUM, KFST4_SUM` | SUM | Radiation time step | main.f90 |
| Geometry areas/volumes | SUM | CC_IBM initialization | geom.f90 |
| CVODE warning cells | SUM | After combustion | fire.f90 |

### 5.5 MPI Communication Flow Through One Time Step

See `graphs/mpi_communication.dot` for the complete graph. Summary:

```
PREDICTOR:
  [Per-mesh computation: INSERT_PARTICLES, VISCOSITY, MASS_TRANSPORT]
  [Per-mesh computation: DENSITY]
      |
      v
  MESH_EXCHANGE(1) ---- species/density exchange ----
      |
  [Per-mesh: VELOCITY_FLUX, WALL_BC, DIVERGENCE_PART_1]
      |
      v
  EXCHANGE_DIVERGENCE_INFO ---- zone pressure integrals (ALLREDUCE) ----
      |
  [Per-mesh: DIVERGENCE_PART_2]
      |
      v
  === PRESSURE_ITERATION_LOOP ===
  |   MESH_EXCHANGE(5) ---- momentum flux + pressure ----
  |   [Per-mesh: PRESSURE_SOLVER]
  |   MESH_EXCHANGE(5) ---- updated pressure ----
  |   [Per-mesh: VELOCITY_ERROR]
  |   MPI_ALLGATHERV ---- velocity/pressure errors ----
  |   (loop until converged)
  ===================================
      |
  [Per-mesh: VELOCITY_PREDICTOR]
      |
      v
  MPI_ALLGATHERV ---- CFL/time step sync ----
      |
  MESH_EXCHANGE(3) ---- velocity/pressure exchange ----
      |
  [Per-mesh: MATCH_VELOCITY, VELOCITY_BC]

CORRECTOR:
  [Per-mesh: VISCOSITY, MASS, DENSITY]
      |
      v
  MESH_EXCHANGE(4) ---- species/density exchange ----
      |
  [Per-mesh: VELOCITY_FLUX, COMBUSTION, PARTICLES]
      |
      v
  MESH_EXCHANGE(7+11) ---- particle exchange ----
      |
  [Per-mesh: WALL_BC, RADIATION]
      |
      v
  MESH_EXCHANGE(2) ---- radiation intensity ----
      |
  [Per-mesh: DIVERGENCE_PART_1]
      |
      v
  EXCHANGE_DIVERGENCE_INFO ---- zone pressure integrals ----
      |
  [Per-mesh: DIVERGENCE_PART_2]
      |
      v
  === PRESSURE_ITERATION_LOOP === (same as predictor)
      |
  [Per-mesh: VELOCITY_CORRECTOR]
      |
      v
  MESH_EXCHANGE(6) ---- final velocity/pressure + wall data ----
      |
  [Per-mesh: MATCH_VELOCITY, VELOCITY_BC, OUTPUT]
      |
      v
  EXCHANGE_GLOBAL_OUTPUTS ---- HRR, mass balance (ALLREDUCE) ----
```

### 5.6 Communication Volume Estimates

Per mesh boundary exchange, the data volume per neighbor:

| Package | Size Formula | Typical Size |
|---------|-------------|--------------|
| Species (CODE 1,4) | NIC * (6 + 2*N_SPECIES) reals | ~thousands of doubles |
| Velocity (CODE 3,6) | IJK_SIZE * 4 reals | ~thousands of doubles |
| Momentum (CODE 5) | NIC * 3 reals + pressure | ~thousands of doubles |
| Radiation (CODE 2) | NRA * N_BANDS * NIC reals | Can be very large |
| Particles (CODE 11) | Variable (per particle) | Depends on count |

Where NIC = number of interface cells between two meshes.

---

## Appendix A: Graph Files

The following DOT graph files are provided in `graphs/`:

| File | Description |
|------|-------------|
| `module_dependencies.dot` | Full module dependency graph |
| `execution_flow.dot` | Time-stepping loop execution order |
| `pressure_iteration.dot` | Pressure solver inner loop detail |
| `data_flow.dot` | Data transformation through one time step |
| `mpi_communication.dot` | Inter-process communication timeline |
| `thread_safety.dot` | Thread safety status of all computation routines |

## Appendix B: POINT_TO_MESH Removal Plan

See [POINT_TO_MESH_REMOVAL.md](POINT_TO_MESH_REMOVAL.md) for a detailed
per-module difficulty assessment and phased removal plan.
