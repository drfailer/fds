# FDS Phase 2: Breaking the Sequential Bottleneck

**Goal:** Reduce sequential fraction from 39% to ~20-25% to achieve 2.0-2.5× overall speedup.

**Current status:** 12 sub-graphs completed, 1.3× overall speedup (9× kernel speedup limited by 39% sequential)

**Phase 2 targets:** CorrRadiation, PredFinal, CorrFinal (potential ~1,000ms savings on 4-mesh test)

---

## Refactoring Pipeline

Each complex routine follows this pipeline:

```
1. ANALYZE ROUTINE
   ├─ Identify OMESH dependencies (cross-mesh reads/writes)
   ├─ Identify local parallelizable loops
   └─ Estimate sequential vs parallel fractions

2. EXTRACT/CONVERT KERNELS
   ├─ Create thread-safe kernels (*_kernels.f90 or inline)
   ├─ Convert callees to TYPE(MESH_TYPE) argument pattern
   └─ Test: compile only

3. CREATE SUB-GRAPH COMPONENTS
   ├─ Data structure (data/*.h)
   ├─ Orchestrator + Collector (state/*.h)
   ├─ Kernel task (task/*.h)
   └─ Dedicated sub-graph wrapper (graph/*_subgraph.h)

4. INTEGRATE INTO MAIN GRAPH
   ├─ Update fds_graph.h
   ├─ Update fds_fortran_interface.h (if needed)
   └─ Test: compile + run test suite

5. VERIFY AND COMMIT
   ├─ All 5 test cases byte-identical
   └─ Commit with descriptive message
```

---

## Phase 2 Roadmap

### Priority 1: VELOCITY_BC Decomposition (Prerequisite for PredFinal/CorrFinal)
**Time:** Embedded in PredFinal (~342ms) and CorrFinal (~385ms)
**Parallelizable:** 75-80% (local wall/edge processing)
**Effort:** High (700+ lines, many edge cases)

**Tasks:**
- [x] Create detailed refactoring plan ([VELOCITY_BC_REFACTORING_PLAN.md](VELOCITY_BC_REFACTORING_PLAN.md))
- [x] **Phase 1**: Analyze CC_VELOCITY callees (SET_GHOSTFACE_VEL_*)
  - [x] Determined module-level variables (FCVAR, CUT_FACE) are mesh-specific
  - [x] **Decision**: Skip Phase 2 - use WallBC pattern with POINT_TO_MESH for module pointers
- [ ] ~~**Phase 2**: Convert CC_VELOCITY callees~~ **SKIPPED**
  - Rationale: WallBC pattern allows POINT_TO_MESH for module pointers (thread-safe per mesh)
- [x] **Phase 3**: Convert VELOCITY_BC to thread-safe VELOCITY_BC_KERNEL ✅
  - [x] Create VELOCITY_BC_KERNEL with `TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M` and `INTEGER, INTENT(IN) :: NM`
  - [x] Call POINT_TO_MESH(NM) at beginning to set up module pointers
  - [x] Use M for explicit accesses (M%WALL, M%U, M%V, etc.)
  - [x] Pass NM to callees when needed
  - [x] Test: full test suite (byte-identical)
  - [x] Commit ✅ (470df20aa7)
- [x] **Phase 4**: Extract VELOCITY_BC components ✅
  - [x] Extract VELOCITY_BC_PREPROCESSING (OMESH wall velocity reads)
  - [x] Extract VELOCITY_BC_PROCESS_EDGES_KERNEL (local edge processing, all edges)
  - [x] Modified VELOCITY_BC_KERNEL to call both components
  - [x] Test: full test suite ✅ (all 5 tests byte-identical)
  - [x] Commit ✅ (638023d8f8)

**Estimated effort**: 11-15 hours total (Phase 2 skipped)
**Actual effort**: 3.5 hours total ✅

**Blockers:** ~~Must complete before PredFinal/CorrFinal refactoring~~ ✅ **UNBLOCKED**

**Current status**: Phase 4 ✅ complete (fds_hh build + tests successful)

**Phase 3 substeps:** ✅ All complete
- [x] Created VELOCITY_BC_CONVERSION_MAP.md with systematic substitution plan
- [x] Created VELOCITY_BC_KERNEL with TYPE(MESH_TYPE) argument
  - Completed ~150 substitutions in 809-line routine
  - Following WallBC pattern: M for explicit accesses, NM for callees, POINT_TO_MESH for module pointers
  - Fixed global arrays: EDGE_COUNT(NM), T_USED(4) remain without M% prefix
- [x] Convert old VELOCITY_BC to wrapper
- [x] Compile test ✅ (successful)
- [x] Full test suite ✅ (all 5 tests byte-identical)
- [x] Commit ✅ (470df20aa7)

**Actual time for Phase 3:** ~1.5 hours (faster than estimated 2-3 hours due to sed automation)

**Test results (Phase 3):**
- ✅ dancing_eddies_1mesh (7.03s) - byte-identical
- ✅ dancing_eddies_2mesh (13.20s) - byte-identical
- ✅ dancing_eddies_4mesh (5.60s) - byte-identical
- ✅ multiple_reac_3mesh (6.14s) - byte-identical
- ✅ species_props_5mesh (1.09s) - byte-identical

**Phase 4 results:** ✅ Complete

Test results (fds_hh):
- ✅ dancing_eddies_1mesh (7.47s) - byte-identical
- ✅ dancing_eddies_2mesh (13.10s) - byte-identical
- ✅ dancing_eddies_4mesh (5.50s) - byte-identical
- ✅ multiple_reac_3mesh (6.11s) - byte-identical
- ✅ species_props_5mesh (1.08s) - byte-identical

Extracted VELOCITY_BC components:
1. **VELOCITY_BC_PREPROCESSING** (~77 lines):
   - WALL_LOOP: reads OMESH velocities for external wall boundaries
   - Initializes M%DRAG_UVWMAX
   - Must run sequentially before parallel kernel (cross-mesh dependencies)

2. **VELOCITY_BC_PROCESS_EDGES_KERNEL** (~762 lines):
   - EDGE_LOOP: processes all cell edges (including INTERPOLATED edges)
   - Parallelizable per-mesh (OMESH reads are safe after preprocessing)
   - INTERPOLATED edges included because M%OMESH data is already synchronized

3. **VELOCITY_BC_KERNEL** (now ~53 lines):
   - Simplified to call preprocessing + process_edges + accounting
   - Maintains backward compatibility
   - Easier to understand and maintain

---

### Priority 2: PredFinal (MATCH_VELOCITY + VELOCITY_BC + SYNTHETIC_TURBULENCE)
**Time:** ~342 ms (4-mesh)
**Parallelizable:** 75-80%
**Effort:** Medium (Fortran ✅ complete, C++ integration pending)

**Status:** Fortran refactoring complete, ready for Hedgehog C++ integration

**Fortran components ready:**
- [x] VELOCITY_BC_PREPROCESSING (sequential - OMESH reads)
- [x] VELOCITY_BC_PROCESS_EDGES_KERNEL (parallel - main computation)
- [ ] MATCH_VELOCITY (already exists - sequential cross-mesh sync)
- [ ] SYNTHETIC_TURBULENCE (exists but needs thread-safe conversion for parallel execution)

**Next: C++ Integration** (see `VELOCITY_BC_HEDGEHOG_INTEGRATION.md`)
- [ ] Option A: Simple task-based integration (2-3 hours)
- [ ] Option B: Dedicated sub-graph with Pattern B (4-6 hours)
- [ ] Test and verify byte-identical
- [ ] Profile performance (target: ~140ms, 2.4× faster)
- [ ] Commit with performance results

**Blockers:** ~~Requires VELOCITY_BC decomposition~~ ✅ **UNBLOCKED** (Fortran work complete)

---

### Priority 3: CorrFinal (MATCH_VELOCITY + VELOCITY_BC + UPDATE_GLOBAL_OUTPUTS)
**Time:** ~385 ms (4-mesh)
**Parallelizable:** 70-75%
**Effort:** Medium (Fortran ✅ complete, C++ integration pending)

**Status:** Fortran refactoring complete, ready for Hedgehog C++ integration

**Fortran components ready:**
- [x] VELOCITY_BC_PREPROCESSING (sequential - OMESH reads)
- [x] VELOCITY_BC_PROCESS_EDGES_KERNEL (parallel - main computation)
- [ ] MATCH_VELOCITY (already exists - sequential cross-mesh sync)
- [ ] UPDATE_GLOBAL_OUTPUTS (exists - needs analysis for parallelization)

**Next: C++ Integration** (see `VELOCITY_BC_HEDGEHOG_INTEGRATION.md`)
- [ ] Integrate with PredFinal or create separate sub-graph
- [ ] Analyze UPDATE_GLOBAL_OUTPUTS for potential parallelization
- [ ] Test and verify byte-identical
- [ ] Profile performance (target: ~170ms, 2.3× faster)
- [ ] Commit with performance results

**Blockers:** ~~Requires VELOCITY_BC decomposition~~ ✅ **UNBLOCKED** (Fortran work complete)

---

---

### Deferred: CorrRadiation (COMPUTE_RADIATION)
**Time:** ~356 ms (4-mesh)
**Parallelizable:** 70-80% (angle loop + spectral bands)
**Effort:** Very High

**Complexity analysis:**
- [x] Analyzed RADIATION_FVM structure
  - Spectral band loops (BAND_LOOP)
  - Angle loops with spatial sweeps (ANGLE_LOOP)
  - Cross-mesh intensity exchange via OMESH%IL_R (INTERPOLATED_BOUNDARY)
  - Global RTE source correction (RAD_Q_SUM, KFST4_SUM accumulation)
  - Multiple geometry cases (cylindrical, 2D, 3D cartesian)

**Why deferred:**
- More complex than initially assessed
- Requires careful handling of angle sweep dependencies
- Cross-mesh intensity interpolation needs sequential barriers
- RTE source correction is global (all meshes contribute)
- Best tackled after gaining experience with VELOCITY_BC decomposition

**Future approach:**
- Parallelize spectral band loop (independent bands)
- Parallelize angle loop within each band (may need careful ordering)
- Keep RTE source correction sequential (global reduction)
- Keep cross-mesh interpolation in preprocessing

**Estimated parallelizable fraction:** 70-80% (revised down from initial 80-85%)

---

## Future Work (Phase 3)

### PressureIteration (Complex Iterative Solver)
**Time:** ~302 ms (2 calls)
**Parallelizable:** 50-60%
**Effort:** Very High (iterative convergence, multiple solver backends)

**Approach:** Hybrid decomposition
- Parallelize per-mesh setup (BAROCLINIC_CORRECTION, PRESSURE_SOLVER_SETUP)
- Keep global solver sequential (FFT/ULMAT/GLMAT)
- Parallelize per-mesh error checking

**Estimated savings:** ~150ms (limited by iterative nature)

---

### Combustion (MPI Load-Balanced)
**Time:** ~100-150 ms
**Parallelizable:** 60-70%
**Effort:** High (MPI coordination)

**Approach:** Per-mesh kernel
- Keep MPI load balancing sequential
- Parallelize COMBUSTION_MODEL calls (expensive ODE integration)
- Requires careful MPI communication pattern

**Estimated savings:** ~60-100ms

---

## Progress Tracking

### Completed (Phase 1)
✅ 12 sub-graphs (VelocityCorrector, VelocityPredictor, DivPart2, CorrStep1, DensityPred, CorrDivPart1, DivSetup, PredStep1, CorrCondens, PredWallDiv, CorrParticle, WallBC)

### Completed (Phase 2)
✅ VELOCITY_BC decomposition (All phases complete - 3.5 hours)
  - Phase 1: CC_VELOCITY analysis ✅
  - Phase 3: Thread-safe conversion ✅ (470df20aa7)
  - Phase 4: Component extraction ✅ (638023d8f8)

### Ready to Start (Phase 2)
🟢 PredFinal (unblocked - VELOCITY_BC components ready)
🟢 CorrFinal (unblocked - VELOCITY_BC components ready)

---

## Performance Targets

### Current (12 sub-graphs, 4-mesh, kernelThreads=4)
- Parallel kernels: 1277 ms (27.2%) — **9.0× speedup achieved**
- Sequential tasks: 1838 ms (39.2%) — **bottleneck**
- Overall: 4.690s → **1.3× speedup vs 1-mesh baseline**

### Target (Phase 2 complete)
- Parallel kernels: ~900 ms (20%) — maintain 9× speedup
- Sequential tasks: ~1100 ms (25%) — **50% reduction**
- Overall: ~3.8s → **2.0-2.5× speedup vs 1-mesh baseline**

**Amdahl's law limit:** With 25% sequential, max speedup = 1/(0.25 + 0.75/9) ≈ **3.3×**

---

## Testing Protocol

After each refactoring step:

1. **Compile:** `cmake --build build_hh --target fds -j$(nproc)`
2. **Run tests:** `cd test_cases && python3 run_tests.py -v`
3. **Verify:** All 5 test cases byte-identical
4. **Commit:** Descriptive message with test results

**Test cases:**
- dancing_eddies_1mesh
- dancing_eddies_2mesh
- multiple_reac_3mesh
- dancing_eddies_4mesh
- species_props_5mesh

---

## Notes

- Follow WallBC pattern: dedicated sub-graph wrapper for traceability
- Document all OMESH access patterns (critical for correctness)
- Maintain byte-identical results at every step
- Ask for help if blocked or uncertain about decomposition strategy
