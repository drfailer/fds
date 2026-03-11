# Hedgehog Velocity Sub-Graph Analysis

## Executive Summary

This document analyzes the current Hedgehog implementation to identify velocity-related tasks and states that would benefit from replacement with lower-level sub-graphs. The goal is to enable **multi-mesh parallel processing within a single node** by having sub-graph tasks process multiple meshes concurrently, while orchestration states manage synchronization and control flow.

## Current High-Level Graph Architecture

### Issues Identified

1. **Sequential mesh processing**: Each task processes one mesh token at a time (numThreads=1 for Phase 1)
2. **Non thread-safe routines**: Many orchestration routines use `POINT_TO_MESH` and global state
3. **Computation in barriers**: Some barrier tasks perform complex computation (e.g., `ChangeTimeStepTask` retry loop)
4. **Fine-grained task granularity**: Many small tasks that could be grouped for parallel execution

### Velocity-Related Tasks

#### Predictor Phase (7 tasks)
1. **PredStep1Task**: `fds_insert_particles`, `fds_compute_viscosity`, `fds_mass_finite_differences`
2. **DensityPredTask**: `fds_density`
3. **PredDivSetupTask**: `fds_set_baroclinic_false`, `fds_viscosity_bc`, `fds_velocity_flux`
4. **PredWallDivTask**: `fds_wall_bc`, `fds_particle_momentum`, `fds_divergence_part_1`
5. **DivPart2PredTask**: `fds_divergence_part_2`
6. **VelPredictorTask**: `fds_velocity_predictor` ← **PRIMARY VELOCITY TASK**
7. **PredFinalTask**: `fds_match_velocity`, `fds_synthetic_turbulence`, `fds_velocity_bc`

#### Corrector Phase (10 tasks)
1. **CorrStep1Task**: `fds_compute_viscosity`, `fds_mass_finite_differences`, `fds_density`
2. **CorrDivSetupTask**: `fds_set_baroclinic_false`, `fds_viscosity_bc`, `fds_velocity_flux`, `fds_agglomeration`
3. **CorrCondensTask**: `fds_condensation`
4. **CorrParticleTask**: `fds_particle_mass_energy`, `fds_move_particles`, `fds_particle_momentum`
5. **CorrWallBCTask**: `fds_wall_bc`
6. **CorrRadiationTask**: `fds_compute_radiation`
7. **CorrDivPart1Task**: `fds_combustion_bc`, `fds_divergence_part_1`
8. **CorrDivPart2Task**: `fds_divergence_part_2`
9. **CorrVelocityTask**: `fds_velocity_corrector`, `fds_check_divergence` ← **PRIMARY VELOCITY TASK**
10. **CorrFinalTask**: `fds_match_velocity`, `fds_velocity_bc`, `fds_update_global_outputs`

## Velocity Routine Analysis

### Thread-Safe Computation Kernels (in velo_kernels.f90)

These routines take `TYPE(MESH_TYPE), INTENT(INOUT) :: M` and avoid `POINT_TO_MESH`:

1. **VELOCITY_FLUX_KERNEL** (velo_kernels.f90:236)
   - Computes convective/diffusive terms for momentum equations
   - Thread-safe: No cross-mesh access, operates on M%
   - Used by: PredDivSetupTask, CorrDivSetupTask

2. **VELOCITY_PREDICTOR_KERNEL** (velo_kernels.f90:121)
   - Updates velocity field: U = U + FVX*DT
   - Thread-safe: No cross-mesh access
   - Used by: VelPredictorTask

3. **VELOCITY_CORRECTOR_KERNEL** (velo_kernels.f90:176)
   - Corrects velocity: U = U + FVX*DT (corrector phase)
   - Thread-safe: No cross-mesh access
   - Used by: CorrVelocityTask

4. **CHECK_STABILITY_KERNEL** (velo_kernels.f90:1315)
   - Computes CFL-limited time step DT_NEW
   - Thread-safe: Writes to DT_NEW(NM), indexed by mesh
   - Used by: CHECK_STABILITY (called after VELOCITY_PREDICTOR)

### Non Thread-Safe Orchestration Routines (in velo.f90)

These routines use `POINT_TO_MESH` and cross-mesh data structures:

1. **VELOCITY_FLUX** (velo.f90:118)
   - **Orchestration**: POINT_TO_MESH, pointer selection (US vs U), CC_IBM handling
   - **Calls kernel**: VELOCITY_FLUX_KERNEL(MESHES(NM), ...)
   - **Thread-safety issues**: POINT_TO_MESH, module-level pointers

2. **VELOCITY_PREDICTOR** (velo.f90:548)
   - **Orchestration**: POINT_TO_MESH, CC_PROJECT_VELOCITY, WALL_VELOCITY_NO_GRADH
   - **Calls kernel**: VELOCITY_PREDICTOR_KERNEL(MESHES(NM), DT)
   - **Calls**: CHECK_STABILITY (writes to global DT_NEW array)
   - **Thread-safety issues**: POINT_TO_MESH, global DT_NEW access pattern

3. **VELOCITY_CORRECTOR** (velo.f90:624)
   - **Orchestration**: POINT_TO_MESH, CC_PROJECT_VELOCITY, WALL_VELOCITY_NO_GRADH
   - **Calls kernel**: VELOCITY_CORRECTOR_KERNEL(MESHES(NM), DT)
   - **Thread-safety issues**: POINT_TO_MESH

4. **VELOCITY_BC** (velo.f90:703)
   - **Heavy cross-mesh access**: EXTERNAL_WALL, OMESH arrays
   - **Not parallelizable**: Must execute after MESH_EXCHANGE
   - **Not a candidate** for sub-graph parallelization

5. **MATCH_VELOCITY** (velo.f90:1518)
   - **Heavy cross-mesh access**: EXTERNAL_WALL, OMESH, MESHES arrays
   - **Not parallelizable**: Interpolates velocities between meshes
   - **Not a candidate** for sub-graph parallelization

## Sub-Graph Candidates (Ranked by Impact)

### 1. VELOCITY_PREDICTOR Sub-Graph (HIGHEST PRIORITY)

**Current Implementation** (fds_graph.h:181-184):
```
predPressureTask → velPredictor → changeTimeStepCollectorSM → changeTimeStepTask → collector3SM
```

**Flow**:
- Sequential: N mesh tokens arrive one at a time to VelPredictorTask
- Each calls `fds_velocity_predictor(t, dt, nm)` which calls `VELOCITY_PREDICTOR_KERNEL`
- All tokens collected by changeTimeStepCollectorSM
- ChangeTimeStepTask checks CFL and may retry entire predictor sequence

**Proposed Sub-Graph Architecture**:

```
[State: VelocityPredictorOrchestrator]
  ├─ Collects all N mesh tokens
  ├─ Emits N parallel work tokens
  │
  ├─> [Task: ParallelVelocityPredictorKernel] (numThreads = N)
  │    ├─ Calls VELOCITY_PREDICTOR_KERNEL(MESHES(NM), DT)
  │    ├─ Calls CHECK_STABILITY_KERNEL(MESHES(NM), DT, DT_NEW, T, NM)
  │    └─ Thread-safe: No POINT_TO_MESH, indexed writes to DT_NEW(NM)
  │
  ├─> [State: CFLCheckOrchestrator]
  │    ├─ Collects all N results
  │    ├─ Checks if ANY(CHANGE_TIME_STEP_INDEX == -1)
  │    ├─ If retry: reduces DT, emits new parallel work tokens → loop back
  │    └─ If success: emits N mesh tokens to next stage
  │
  └─> [Output: N mesh tokens with updated velocities and DT_NEW]
```

**Benefits**:
- **Parallel execution**: All N meshes process velocity simultaneously (N threads)
- **Eliminates POINT_TO_MESH**: Kernel takes MESHES(NM) directly
- **Handles CFL retry**: Orchestrator state manages retry loop without re-entering graph
- **Major speedup potential**: Currently sequential bottleneck

**Challenges**:
- **CC_PROJECT_VELOCITY**: Called in orchestration layer (velo.f90:576), may not be thread-safe
- **WALL_VELOCITY_NO_GRADH**: Called for ULMAT pressure solvers (velo.f90:580), thread-safety unknown
- **DT_NEW array**: Indexed writes are safe, but retry logic accesses full array

**Mitigation**:
- Move CC_IBM special handling to sequential orchestrator state (before or after parallel kernel)
- For WALL_VELOCITY_NO_GRADH, assess thread-safety or serialize this step
- CFLCheckOrchestrator handles global DT_NEW reduction sequentially

---

### 2. VELOCITY_CORRECTOR Sub-Graph (HIGH PRIORITY)

**Current Implementation** (fds_graph.h:218-220):
```
corrPressureTask → corrVelocity → collector6bSM → meshExchange6b
```

**Flow**:
- Sequential: N mesh tokens arrive one at a time to CorrVelocityTask
- Each calls `fds_velocity_corrector(t, dt, nm)` and `fds_check_divergence(nm)`
- All tokens collected for MESH_EXCHANGE(6)

**Proposed Sub-Graph Architecture**:

```
[State: VelocityCorrectorOrchestrator]
  ├─ Collects all N mesh tokens
  ├─ Emits N parallel work tokens
  │
  ├─> [Task: ParallelVelocityCorrectorKernel] (numThreads = N)
  │    ├─ Calls VELOCITY_CORRECTOR_KERNEL(MESHES(NM), DT)
  │    ├─ Calls CHECK_DIVERGENCE_KERNEL(MESHES(NM)) [if extracted to kernel]
  │    └─ Thread-safe: No POINT_TO_MESH
  │
  ├─> [State: DivergenceCheckCollector]
  │    ├─ Collects all N results
  │    └─ Emits N mesh tokens to next stage
  │
  └─> [Output: N mesh tokens with corrected velocities]
```

**Benefits**:
- **Parallel execution**: All N meshes process velocity correction simultaneously
- **Eliminates POINT_TO_MESH**: Kernel takes MESHES(NM) directly
- **Simpler than predictor**: No CFL retry loop to handle

**Challenges**:
- **CC_PROJECT_VELOCITY**: Called before and after kernel (velo.f90:648, 662), may not be thread-safe
- **WALL_VELOCITY_NO_GRADH**: Called for ULMAT solvers (velo.f90:653, 667)
- **CHECK_DIVERGENCE**: Currently not extracted to kernel, may use POINT_TO_MESH

**Mitigation**:
- Move CC_IBM special handling to sequential orchestrator state
- Extract CHECK_DIVERGENCE to kernel (similar to existing divergence kernels)
- Serialize WALL_VELOCITY_NO_GRADH or assess thread-safety

---

### 3. VELOCITY_FLUX Sub-Graph (MEDIUM PRIORITY)

**Current Implementation**:
- Called within PredDivSetupTask (fds_graph.h:38, predictor_tasks.h:51)
- Called within CorrDivSetupTask (fds_graph.h:46, corrector_tasks.h:35)

**Flow**:
- Sequential: Each task calls `fds_velocity_flux(t, dt, nm, estimated)`
- Kernel VELOCITY_FLUX_KERNEL is already thread-safe

**Proposed Sub-Graph Architecture**:

```
[State: VelocityFluxOrchestrator]
  ├─ Collects all N mesh tokens
  ├─ Emits N parallel work tokens
  │
  ├─> [Task: ParallelVelocityFluxKernel] (numThreads = N)
  │    ├─ Calls VELOCITY_FLUX_KERNEL(MESHES(NM), T, DT, NM, APPLY_TO_ESTIMATED_VARIABLES, GX, GY, GZ)
  │    └─ Thread-safe: No POINT_TO_MESH in kernel
  │
  └─> [Output: N mesh tokens with updated velocity fluxes (FVX, FVY, FVZ)]
```

**Benefits**:
- **Parallel execution**: Velocity flux computation for all meshes simultaneously
- **Kernel already thread-safe**: Main computation is in VELOCITY_FLUX_KERNEL

**Challenges**:
- **CC_IBM handling**: CUTFACE_VELOCITIES and CC_VELOCITY_FLUX calls (velo.f90:154, 170)
- **VELOCITY_FLUX_CYLINDRICAL**: Separate path for cylindrical coordinates (velo.f90:186)
- **Embedded in larger tasks**: Currently bundled with viscosity_bc in DivSetup tasks

**Mitigation**:
- Move CC_IBM special handling to sequential orchestrator
- Keep cylindrical path separate or create cylindrical sub-graph variant
- Extract velocity_flux into dedicated sub-graph (unbundle from DivSetup)

---

## Recommended Implementation Strategy

### Phase 1: Single Sub-Graph Prototype

**Target**: VELOCITY_CORRECTOR sub-graph
- **Rationale**: Simpler than predictor (no retry loop), significant parallel benefit
- **Steps**:
  1. Create VelocityCorrectorOrchestrator state (collects N tokens)
  2. Create ParallelVelocityCorrectorKernel task (numThreads = N)
  3. Extract CHECK_DIVERGENCE to kernel (if not already done)
  4. Handle CC_IBM/WALL_VELOCITY_NO_GRADH in orchestrator (sequential pre/post steps)
  5. Replace CorrVelocityTask + collector with sub-graph
  6. Test with numThreads=1 (sequential, byte-identical results)
  7. Test with numThreads=N (parallel, byte-identical results if kernels truly thread-safe)

### Phase 2: Predictor Sub-Graph

**Target**: VELOCITY_PREDICTOR sub-graph with CFL retry
- **Rationale**: Highest impact on performance, more complex due to retry logic
- **Steps**:
  1. Create VelocityPredictorOrchestrator state
  2. Create ParallelVelocityPredictorKernel task (includes CHECK_STABILITY_KERNEL)
  3. Create CFLCheckOrchestrator state (manages retry loop)
  4. Replace VelPredictorTask + ChangeTimeStepTask with sub-graph
  5. Test retry logic with intentionally small DT to trigger retries
  6. Validate byte-identical results

### Phase 3: Additional Sub-Graphs

**Targets**: VELOCITY_FLUX, other computational kernels (DIVERGENCE, MASS, etc.)
- **Rationale**: Apply same pattern to other kernel-based operations
- **Steps**: Similar to Phase 1 for each kernel module

---

## Thread-Safety Requirements Checklist

For each sub-graph to enable parallel multi-mesh processing:

### Must Eliminate:
- ✗ **POINT_TO_MESH calls**: Replaced by passing MESHES(NM) to kernels
- ✗ **Module-level pointer aliases**: Use explicit M%ARRAY indexing in kernels
- ✗ **Unindexed global writes**: DT_NEW must be DT_NEW(NM)

### Must Handle Carefully:
- ⚠ **Indexed global arrays**: DT_NEW(NM), CHANGE_TIME_STEP_INDEX(NM) - safe if each thread writes unique index
- ⚠ **Read-only global data**: GLOBAL_CONSTANTS - safe
- ⚠ **CC_IBM special handling**: Likely not thread-safe, move to sequential orchestrator
- ⚠ **Cross-mesh operations**: MESH_EXCHANGE, OMESH - must remain in barriers between sub-graphs

### Must Verify:
- ❓ **WALL_VELOCITY_NO_GRADH**: Thread-safety unknown, investigate pois.f90
- ❓ **CC_PROJECT_VELOCITY**: Thread-safety unknown, investigate ccib_velocity.f90
- ❓ **CHECK_DIVERGENCE**: Thread-safety unknown, may need kernel extraction

---

## Expected Performance Impact

### Current Phase 1 (Sequential, numThreads=1):
- Each mesh token processed one at a time through velocity tasks
- Total velocity time per time step: `N_meshes × T_velocity_per_mesh`

### Proposed Phase 2 (Parallel, numThreads=N):
- All mesh tokens processed simultaneously in sub-graph tasks
- Total velocity time per time step: `T_velocity_per_mesh + T_orchestration_overhead`
- **Theoretical speedup**: ~N× for velocity operations (significant portion of time step)

### Realistic Speedup Estimate:
- If velocity operations are 30% of time step, N=4 meshes: ~2.4× overall speedup
- If velocity operations are 50% of time step, N=4 meshes: ~3.2× overall speedup
- Diminishing returns beyond N=number of CPU cores

---

## Implementation Notes

### Sub-Graph Pattern

All velocity sub-graphs follow this pattern:

1. **Orchestrator State** (pure data-flow):
   - Collects N MeshData tokens
   - Performs any sequential pre-processing (CC_IBM setup, etc.)
   - Emits N work tokens to parallel task

2. **Parallel Kernel Task** (computation, numThreads=N):
   - Receives work token (contains mesh index NM)
   - Calls thread-safe kernel: `KERNEL(MESHES(NM), ...)`
   - No POINT_TO_MESH, no cross-mesh access
   - Emits result token

3. **Collector State** (pure data-flow):
   - Collects N result tokens
   - Performs any sequential post-processing (diagnostics, CFL check, etc.)
   - Decides: continue (emit N MeshData) or retry (emit new work tokens)

### C++ Hedgehog Implementation

```cpp
// Orchestrator state (collects and dispatches)
class VelocityPredictorOrchestrator : public hh::AbstractState<1, MeshData, VelWorkToken> {
public:
    VelocityPredictorOrchestrator(int nmeshes) : nmeshes_(nmeshes) {
        collected_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<MeshData> data) override {
        collected_.push_back(data);
        if (collected_.size() == nmeshes_) {
            // Sequential pre-processing here (CC_IBM, etc.)
            for (auto &md : collected_) {
                auto work = std::make_shared<VelWorkToken>(md->nm, md->t, md->dt);
                this->addResult(work);
            }
            collected_.clear();
        }
    }
private:
    int nmeshes_;
    std::vector<std::shared_ptr<MeshData>> collected_;
};

// Parallel kernel task
class ParallelVelocityPredictorKernel : public hh::AbstractTask<1, VelWorkToken, VelResultToken> {
public:
    ParallelVelocityPredictorKernel(size_t numThreads)
        : hh::AbstractTask<1, VelWorkToken, VelResultToken>("ParallelVelPred", numThreads) {}

    void execute(std::shared_ptr<VelWorkToken> work) override {
        // Thread-safe kernel call
        fds_velocity_predictor_kernel(work->nm, work->t, work->dt);
        fds_check_stability_kernel(work->nm, work->dt);

        auto result = std::make_shared<VelResultToken>(work->nm);
        this->addResult(result);
    }

    std::shared_ptr<hh::AbstractTask<1, VelWorkToken, VelResultToken>> copy() override {
        return std::make_shared<ParallelVelocityPredictorKernel>(this->numberThreads());
    }
};

// CFL check collector state
class CFLCheckOrchestrator : public hh::AbstractState<1, VelResultToken, MeshData> {
public:
    CFLCheckOrchestrator(int nmeshes) : nmeshes_(nmeshes) {
        results_.reserve(nmeshes);
    }

    void execute(std::shared_ptr<VelResultToken> result) override {
        results_.push_back(result);
        if (results_.size() == nmeshes_) {
            // Check CFL globally (sequential)
            int needRetry = 0;
            double newDt = 0.0;
            fds_check_change_time_step(&needRetry, &newDt);

            if (needRetry) {
                // Emit new work tokens with reduced DT → loop back to parallel task
                // (Complex retry logic from ChangeTimeStepTask)
            } else {
                // Success: emit mesh tokens to next stage
                for (auto &md : originalMeshData_) {
                    this->addResult(md);
                }
            }
            results_.clear();
        }
    }
private:
    int nmeshes_;
    std::vector<std::shared_ptr<VelResultToken>> results_;
    std::vector<std::shared_ptr<MeshData>> originalMeshData_;
};
```

---

## Conclusion

The **VELOCITY_PREDICTOR** and **VELOCITY_CORRECTOR** sub-graphs represent the highest-value targets for enabling multi-mesh parallel processing within the Hedgehog framework. These operations:

1. Are currently sequential bottlenecks (process one mesh at a time)
2. Have thread-safe computation kernels already extracted
3. Consume significant wall-clock time per time step
4. Can be parallelized with orchestrator states managing synchronization

The recommended approach is to prototype the VELOCITY_CORRECTOR sub-graph first (simpler, no retry logic), validate byte-identical results and performance gains, then apply the same pattern to VELOCITY_PREDICTOR (with CFL retry orchestration) and eventually other kernel-based modules.

The sub-graph pattern (Orchestrator State → Parallel Kernel Task → Collector State) provides a clean separation between:
- **Data-flow control** (states)
- **Parallel computation** (tasks calling thread-safe kernels)
- **Cross-mesh synchronization** (barriers between sub-graphs)

This architecture aligns with Hedgehog's strengths and FDS's kernel extraction work, enabling significant performance improvements for multi-mesh simulations on modern multi-core nodes.
