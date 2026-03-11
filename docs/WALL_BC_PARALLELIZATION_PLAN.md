# WALL_BC Parallelization Implementation Plan

## Status: Ready to Implement

All prerequisite thread-safe conversions completed:
- ✅ CALC_HVAC_BC (52 lines)
- ✅ HEAT_TRANSFER_COEFFICIENT (~175 lines)
- ✅ SURFACE_HEAT_TRANSFER (379 lines)
- ✅ CALCULATE_ZZ_F (413 lines)

## Three-Phase Architecture

### Phase 1: ASSIGN_GHOST_VALUE (Sequential)
**Location**: Already extracted in wall.f90
**Why sequential**: Reads from OMESH (neighboring mesh data)
**Cells processed**: ~5-10% (INTERPOLATED_BOUNDARY cells only)

**Current implementation**:
```fortran
SUBROUTINE ASSIGN_GHOST_VALUE(NM,T,DT)
  ! Loop over EXTERNAL_WALL cells with INTERPOLATED_BOUNDARY
  ! Reads neighboring mesh via OMESH (cross-mesh access)
  ! Sets ghost cell values
END SUBROUTINE
```

### Phase 2: WALL_BC_PROCESS_CELLS_KERNEL (Parallel)
**To be created**: New kernel in wall_kernels.f90
**Why parallelizable**: Cell-local processing, no cross-mesh access
**Cells processed**: ~90% (all cells WITHOUT HAS_INTERPOLATED_BC or HAS_BACK_MESH flags)

**Proposed signature**:
```fortran
SUBROUTINE WALL_BC_PROCESS_CELLS_KERNEL(M, PREDICTOR_FLAG, T, DT, CALL_HT_1D)
  TYPE(MESH_TYPE), POINTER :: M
  LOGICAL, INTENT(IN) :: PREDICTOR_FLAG, CALL_HT_1D
  REAL(EB), INTENT(IN) :: T, DT

  ! Loop over wall cells (IW=1 to N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS)
  ! Skip cells with WC%HAS_INTERPOLATED_BC or WC%HAS_BACK_MESH

  ! Calls thread-safe routines:
  ! - SURFACE_HEAT_TRANSFER(NM, PREDICTOR_FLAG, ...)
  ! - CALCULATE_ZZ_F(NM, PREDICTOR_FLAG, ...)
  ! - CALCULATE_RHO_F_KERNEL(M, ...)
  ! - NEAR_SURFACE_GAS_VARIABLES_KERNEL(M, ...)

  ! Similar loop for CFACE cells
  ! Similar loop for particles
END SUBROUTINE
```

### Phase 3: Cross-Mesh Finalization (Sequential)
**To be created**: New orchestration routine
**Why sequential**: Handles BACK_MESH coupling and CONSUME_MASS (updates neighboring mesh)
**Cells processed**: ~5% (thin walls spanning meshes, particle off-gassing)

**Proposed implementation**:
```fortran
SUBROUTINE WALL_BC_FINALIZE(NM,T,DT)
  ! Process cells with HAS_BACK_MESH flag
  !   - SOLID_HEAT_TRANSFER with BACK_MESH > 0
  !   - Uses MESHES(BACK_MESH) access

  ! Process particle off-gassing
  !   - DEPOSIT_PARTICLE_MASS (calls fds_deposit_mass_to_mesh)
  !   - Updates neighboring mesh via OMESH
END SUBROUTINE
```

## Implementation Steps

### Step 1: Extract WALL_BC_PROCESS_CELLS_KERNEL

1. **Create kernel in wall_kernels.f90**:
   ```fortran
   SUBROUTINE WALL_BC_PROCESS_CELLS_KERNEL(M, PREDICTOR_FLAG, T, DT, CALL_HT_1D)
   ```

2. **Extract from wall.f90 lines ~140-270**:
   - Wall cell loop (skip cells with HAS_INTERPOLATED_BC or HAS_BACK_MESH)
   - CFACE loop (similar filtering)
   - Particle loop (skip cells requiring DEPOSIT_PARTICLE_MASS)

3. **Add module-level pointers** (like SURFACE_HEAT_TRANSFER pattern):
   ```fortran
   TYPE(MESH_TYPE), POINTER :: M
   REAL(EB), POINTER, DIMENSION(:,:,:) :: UU, VV, WW, RHOP
   REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP

   M => MESHES(NM)
   IF (PREDICTOR_FLAG) THEN
      UU => M%US; VV => M%VS; WW => M%WS; RHOP => M%RHOS; ZZP => M%ZZS
   ELSE
      UU => M%U; VV => M%V; WW => M%W; RHOP => M%RHO; ZZP => M%ZZ
   ENDIF
   ```

4. **Replace module pointer accesses** with M% prefix

5. **Build and test**: Verify byte-identical results

### Step 2: Create WALL_BC_FINALIZE Orchestration

1. **Extract finalization code** from WALL_BC:
   - BACK_MESH handling from SOLID_HEAT_TRANSFER calls
   - DEPOSIT_PARTICLE_MASS calls

2. **Keep sequential** (no kernel extraction needed)

### Step 3: Restructure Main WALL_BC

Update wall.f90 WALL_BC to three-phase structure:

```fortran
SUBROUTINE WALL_BC(T,DT,NM)
  ! Phase 1: Sequential cross-mesh reads
  IF (N_EXTERNAL_WALL_CELLS > 0) CALL ASSIGN_GHOST_VALUE(NM,T,DT)

  ! Phase 2: Parallel cell processing (via hedgehog)
  CALL WALL_BC_PROCESS_CELLS_KERNEL(MESHES(NM), PREDICTOR, T, DT, CALL_HT_1D)

  ! Phase 3: Sequential cross-mesh writes
  CALL WALL_BC_FINALIZE(NM,T,DT)
END SUBROUTINE
```

### Step 4: Create Hedgehog Sub-Graph

**File**: `Source/hedgehog/graph/wallbc_subgraph.h`

```cpp
auto wallBCGraph = std::make_shared<hh::Graph<...>>("WallBC");

// Orchestrator (Phase 1): run ASSIGN_GHOST_VALUE sequentially
auto wallBCOrchestrator = std::make_shared<WallBCOrchestratorState>(...);

// Parallel kernel (Phase 2): run WALL_BC_PROCESS_CELLS_KERNEL
auto wallBCKernelTask = std::make_shared<WallBCKernelTask>(kernelThreads);

// Collector + Finalization (Phase 3): run WALL_BC_FINALIZE sequentially
auto wallBCCollector = std::make_shared<WallBCCollectorState>(...);

wallBCGraph->input<MeshData>(wallBCOrchestrator);
wallBCGraph->edges(
    wallBCOrchestrator, wallBCKernelTask,
    wallBCKernelTask, wallBCCollector
);
wallBCGraph->output<MeshData>(wallBCCollector);
```

### Step 5: Create C Wrappers

**File**: `Source/fds_c_interface.f90`

```fortran
RECURSIVE SUBROUTINE fds_wall_bc_process_cells_kernel(nm, predictor_flag, t, dt, call_ht_1d) &
    BIND(C, NAME='fds_wall_bc_process_cells_kernel')
  INTEGER(C_INT), INTENT(IN), VALUE :: nm
  LOGICAL(C_BOOL), INTENT(IN), VALUE :: predictor_flag, call_ht_1d
  REAL(C_DOUBLE), INTENT(IN), VALUE :: t, dt
  CALL WALL_BC_PROCESS_CELLS_KERNEL(MESHES(nm), predictor_flag, t, dt, call_ht_1d)
END SUBROUTINE
```

**File**: `Source/hedgehog/fds_fortran_interface.h`

```cpp
extern "C" {
  void fds_wall_bc_process_cells_kernel(int nm, bool predictor_flag,
                                         double t, double dt, bool call_ht_1d);
}
```

### Step 6: Integration Testing

1. **Build both targets**: `fds` and `fds_hh`
2. **Run 1-mesh test**: Verify byte-identical results
3. **Run 4-mesh test**: Verify byte-identical DEVC (HRR may have FP noise)
4. **Performance check**: Compare wall clock time (expect speedup with multiple meshes)

## Expected Outcomes

**Parallelization gain**: ~90% of WALL_BC processing can run concurrently across meshes

**Thread safety**: All callees already converted (SURFACE_HEAT_TRANSFER, CALCULATE_ZZ_F, etc.)

**Maintainability**: Clear three-phase structure maps directly to graph nodes

**Next targets**: After WALL_BC, tackle RADIATION and other sequential bottlenecks

## Files to Modify

1. `Source/wall_kernels.f90` - Add WALL_BC_PROCESS_CELLS_KERNEL
2. `Source/wall.f90` - Add WALL_BC_FINALIZE, restructure WALL_BC
3. `Source/fds_c_interface.f90` - Add C wrapper
4. `Source/hedgehog/fds_fortran_interface.h` - Add C declaration
5. `Source/hedgehog/graph/wallbc_subgraph.h` - Create sub-graph (NEW FILE)
6. `Source/hedgehog/data/wallbc_data.h` - Create data structures (NEW FILE)
7. `Source/hedgehog/state/wallbc_state.h` - Create states (NEW FILE)
8. `Source/hedgehog/task/wallbc_kernel_task.h` - Create kernel task (NEW FILE)
9. `Source/hedgehog/graph/fds_graph.cpp` - Integrate WallBC sub-graph

## Estimated Effort

**Kernel extraction**: 2-3 hours
**Hedgehog integration**: 2-3 hours
**Testing and debugging**: 1-2 hours

**Total**: 5-8 hours to fully parallelized WALL_BC
