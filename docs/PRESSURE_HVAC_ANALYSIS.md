# Pressure Iteration & HVAC Parallelization Analysis

## HVAC_CALC — No parallelization possible

**Location**: `Source/hvac.f90:1370-1523`
**C wrapper**: `fds_hvac_calc` in `fds_c_interface.f90:636-640`

HVAC_CALC is a pure global network solver:
- **Zero per-mesh loops** — all mesh data is pre-aggregated into global `NODE_PROPERTIES` arrays by `HVAC_BC_IN(NM)` (called separately per-mesh before HVAC_CALC)
- Runs on **RANK==0 only** in MPI mode
- Solves a coupled system of duct/node equations via Newton-Raphson iteration with Gaussian elimination
- N_NETWORKS are looped sequentially (typically 1-10 networks)

**Per-mesh work**: Only `HVAC_BC_IN(NM)` is per-mesh, and it's already called in a separate DO NM loop before HVAC_CALC.

**Conclusion**: Must remain a barrier task. No kernel extraction or sub-graph possible.

## PRESSURE_ITERATION_SCHEME — Per-mesh loops inside convergence loop

**Location**: `Source/hedgehog/fds_driver.f90:568-693`
**C wrapper**: `fds_pressure_iteration` in `fds_c_interface.f90:585-588`

### Structure per iteration

```
PRESSURE_ITERATION_LOOP: DO              ← sequential (convergence-dependent)

  [conditional: iteration 1 or ITERATE_BAROCLINIC_TERM]
    Per-mesh: BAROCLINIC_CORRECTION        ✅ parallelizable
    Per-mesh: CC_NO_FLUX(.TRUE.)           ✅ parallelizable
    MESH_EXCHANGE(5)                       ❌ global barrier
    Per-mesh: MATCH_VELOCITY_FLUX          ⚠️ reads OMESH (after exchange)

  Per-mesh: NO_FLUX                        ⚠️ reads OMESH (after exchange)
  Per-mesh: CC_NO_FLUX(.FALSE.)            ✅ parallelizable
  Per-mesh: PRESSURE_SOLVER_COMPUTE_RHS    ✅ parallelizable

  SELECT CASE(PRES_FLAG)
    FFT:   Per-mesh PRESSURE_SOLVER_FFT    ✅ parallelizable
    ULMAT: Per-mesh ULMAT_SOLVER           ✅ parallelizable
    GLMAT: GLMAT_SOLVER (global)           ❌ coupled multi-mesh matrix
           MESH_EXCHANGE(5)               ❌ global barrier
           COPY_H_OMESH_TO_MESH           ❌ global

  Per-mesh: CHECK_RESIDUALS                ✅ parallelizable

  IF (.NOT.ITERATE_PRESSURE) EXIT

  MESH_EXCHANGE(5)                         ❌ global barrier
  Per-mesh: COMPUTE_VELOCITY_ERROR         ⚠️ reads OMESH (after exchange)
  Per-mesh: CC_COMPUTE_VELOCITY_ERROR      ⚠️ reads OMESH (after exchange)

  MPI_ALLGATHERV                           ❌ global reduction
  Convergence checks (MAXVAL)              ❌ global
  → EXIT or continue

ENDDO PRESSURE_ITERATION_LOOP
```

### Per-mesh routines

| Subroutine | Cross-mesh? | Notes |
|------------|-------------|-------|
| BAROCLINIC_CORRECTION(T,NM) | No | Pure per-mesh |
| CC_NO_FLUX(DT,NM,FLAG) | No | Pure per-mesh |
| MATCH_VELOCITY_FLUX(NM) | Yes (OMESH FVX/FVY/FVZ) | After MESH_EXCHANGE(5) |
| NO_FLUX(DT,NM) | Yes (OMESH H/HS) | After MESH_EXCHANGE(5) |
| PRESSURE_SOLVER_COMPUTE_RHS(T,DT,NM) | No | Pure per-mesh |
| PRESSURE_SOLVER_FFT(NM) | No | Pure per-mesh |
| ULMAT_SOLVER(NM,T,DT) | No | Pure per-mesh |
| CHECK_RESIDUALS(NM) | No | Writes PRESSURE_ERROR_MAX(NM) |
| COMPUTE_VELOCITY_ERROR(DT,NM) | Yes (OMESH U/V/W/H) | After MESH_EXCHANGE(5) |

### Global (non-parallelizable) routines

| Subroutine | Reason |
|------------|--------|
| GLMAT_SOLVER(T,DT) | Coupled multi-mesh pressure system |
| TUNNEL_POISSON_SOLVER | Global preconditioner |
| MESH_EXCHANGE(5) | Global synchronization (2-3 per iteration) |
| MPI_ALLGATHERV | Global reduction of error arrays |
| COPY_H_OMESH_TO_MESH | Distributes GLMAT solution |

### Constraints

1. Outer loop **cannot be unrolled** — exit depends on global convergence (MAXVAL across all meshes)
2. **2-3 MESH_EXCHANGE(5)** barriers per iteration
3. GLMAT case is inherently global (coupled matrix across meshes)
4. Multiple conditional paths (CC_IBM, solver type, baroclinic)
5. Variable iteration count (data-dependent convergence)

### Profiling

From 4-mesh dancing_eddies (27 timesteps):
- PressureIteration (2× per timestep): **302ms total** (~5.6ms per call)
- Per-mesh work estimated at 60-70% of per-iteration compute

### Possible approach (deferred)

A pressure iteration sub-graph with an internal cycle (similar to changeTimeStepSubgraph) could parallelize the per-mesh loops within each iteration:

```
Entry(BarrierData)
  → BaroclinicKernel (parallel, conditional)
  → MeshExchange5 (barrier)
  → NoFlux+RHS Kernel (parallel)
  → PressureSolver Kernel (parallel, solver-dependent)
  → CheckResiduals Kernel (parallel)
  → [if ITERATE_PRESSURE]
      → MeshExchange5 (barrier)
      → VelocityError Kernel (parallel)
      → ConvergenceCheck → cycle or exit
Exit(MeshData)
```

Complexity is high due to conditional paths and multiple internal barriers. Benefit is modest at current scale (~5.6ms per call). May become worthwhile with larger meshes or more meshes.
