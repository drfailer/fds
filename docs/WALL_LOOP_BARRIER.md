# Wall Loop Barrier to Intra-Mesh Parallelism

## Problem Statement

Block decomposition along K partitions a mesh into sub-blocks `[K1:K2]` that can be processed
in parallel by multiple threads. This works for pure I,J,K cell loops where each cell's update
is independent. However, **13 out of 16 kernel tasks** cannot use this approach because they
contain **wall loops** — iterations over wall cell indices `IW=1 to N_WALL_CELLS` — that
access arbitrary `(I,J,K)` cells determined at runtime by geometry.

### Why Wall Loops Block K-Decomposition

Wall cells are indexed by a flat integer `IW`, not by `(I,J,K)`. Each wall cell `IW` maps to
a gas cell at `(IIG, JJG, KKG)` via `BOUNDARY_COORD(WC%BC_INDEX)`. A wall at `K=5` and a wall
at `K=35` can both appear in the same wall loop. This means:

1. **Cannot partition by K**: A K-block `[K1:K2]` doesn't know which wall cells belong to it
   without pre-filtering the wall list. Wall cells at any K can read/write arrays at any K.

2. **Race conditions**: If two threads process different K-blocks but both contain wall cells
   that write to the same array elements (e.g., adjacent gas cells of a thin wall), data races
   occur.

3. **Sequential dependencies**: Many wall loops modify the same interim arrays that the
   preceding cell loop computed (e.g., `CONDENSATION_EVAPORATION_KERNEL` alternates
   cell loops and wall loops within a species iteration, sharing `ZZ_INTERIM`).

### Performance Impact

Benchmarking shows the wall-loop-containing kernels dominate runtime:

```
Kernel                    % of Runtime    Blocking Issue
─────────────────────────────────────────────────────────
PredWallDivKernel            11.8%        DIVERGENCE_PART_1 has wall loops + accumulators
CorrDivPart1Kernel            9.5%        Same kernel as above (corrector instance)
WallBCKernel (pred)           8.8%        Pure wall loop (IW=1 to N_WALL_CELLS)
WallBCKernel (corr)           7.4%        Same
CorrRadiationKernel           8.2%        Wall loops + particle loops + angle sweeps
VelocityBCEdges (pred)        7.3%        External wall loop + edge loop
VelocityBCEdges (corr)        7.2%        Same
CorrStep1Kernel               6.9%        COMPUTE_VISCOSITY wall loops + MASS_FINITE_DIFFS wall loops
DivSetupKernel (pred)         6.1%        VELOCITY_FLUX wall loops
DivSetupKernel (corr)         5.0%        Same
DumpMeshOutputs               4.8%        I/O (inherently sequential)
PredStep1Kernel               4.3%        Same wall loops as CorrStep1
PressureIteration (pred)      2.7%        FFT solver (global mesh operation)
PressureIteration (corr)      2.1%        Same
DensityPredKernel             1.2%        Zone loops + wall loops + CHECK_MASS_DENSITY
DivPart2Kernel (×2)           2.3%        Zone loops + wall loops
─────────────────────────────────────────────────────────
SUBTOTAL (blocked)           95.5%        Cannot use K-block decomposition

VelPredBlockKernel            0.1%        ✓ Block-decomposed
VelCorrBlockKernel            0.1%        ✓ Block-decomposed
PartMomBlockKernel            0.0%        ✓ Block-decomposed
CheckStability/CheckDiv       0.6%        Mesh-level reductions (post-block)
Other (barriers, I/O, etc.)   3.7%        Sequential infrastructure
```

The 3 block-decomposed kernels account for **0.2% of runtime** — Amdahl's law limits
the speedup to essentially zero.

## Detailed Analysis: Wall Loops Per Kernel

### 1. DIVERGENCE_PART_1_KERNEL (divg_kernels.f90:27-1387) — 21.3% of runtime

Used by: PredWallDivKernelTask, CorrDivPart1KernelTask, RetryMomentumDivKernelTask

**Wall loops:**
- `WALL_LOOP3` (line 99): `IW=1 to N_EXTERNAL+N_INTERNAL` — zeros `Q_LEAK` on `OPEN_BOUNDARY` walls
- `WALL_LOOP` (line 195): `IW=1 to N_EXTERNAL+N_INTERNAL` — computes wall enthalpy flux contributions to divergence. Reads `B1%Q_RAD_IN/OUT`, writes `M%D_SOURCE(IIG,JJG,KKG)`, `M%M_DOT_PPP(IIG,JJG,KKG,:)`. **This is the critical one** — it writes to gas cells at arbitrary (I,J,K) locations.
- `WALL_LOOP_2` (line 315): `IW=1 to N_EXTERNAL+N_INTERNAL` — computes wall species mass flux. Writes `M%D_SOURCE`, `M%M_DOT_PPP` at gas cells.
- `BOUNDARY_LOOP` (line 469): `IW=1 to N_EXTERNAL` — applies OPEN/MIRROR/INTERPOLATED boundary corrections
- `CORRECTION_LOOP` (line 515): `IW=1 to N_EXTERNAL+N_INTERNAL` — density correction at boundary
- `WALL_LOOP4` (line 734): `IW=1 to N_EXTERNAL+N_INTERNAL` — pressure zone wall contributions to `D_SUM_LOC`/`P_SUM_LOC`
- `WALL_LOOP` (line 807): `IW=1 to N_EXTERNAL+N_INTERNAL` — `DPVOL` wall contribution
- `WALL_LOOP_2` (line 988): `IW=1 to N_EXTERNAL+N_INTERNAL` — mass flux wall contribution
- `WALL_LOOP_3` (line 1081): `IW=1 to N_EXTERNAL+N_INTERNAL` — species diffusion wall contribution

**Cell loops between wall loops:** Multiple I,J,K loops interspersed with the wall loops above, sharing the same arrays (`D_SOURCE`, `M_DOT_PPP`, divergence work arrays).

**Key difficulty:** The wall loops are not at the beginning or end — they are **interleaved** throughout the kernel, modifying the same arrays that cell loops read/write. Cannot simply split into "wall pre-processing → cell kernel → wall post-processing."

### 2. WALL_BC_PROCESS_CELLS_KERNEL — 16.2% of runtime

Used by: WallBCKernelTask (predictor and corrector instances)

**Structure:** Single wall loop `IW=1 to N_WALL_CELLS` calling `SURFACE_HEAT_TRANSFER`, `CALCULATE_ZZ_F`, etc. per wall cell. This is inherently a wall-indexed operation.

### 3. VELOCITY_BC_PROCESS_EDGES_KERNEL (velo_kernels.f90:1430) — 14.5% of runtime

Used by: VelocityBCEdgesTask (predictor and corrector instances)

**Wall loops:**
- `WALL_LOOP` (line 1425 in VISCOSITY_BC_KERNEL, called before this): `IW=1 to N_EXTERNAL+N_INTERNAL`
- The edges kernel itself loops over mesh edges applying wall BC conditions

### 4. COMPUTE_VISCOSITY_KERNEL (velo_kernels.f90:921) — 11.2% of runtime

Used by: PredStep1KernelTask, CorrStep1KernelTask

**Wall loops:**
- `DO IW=1,N_EXTERNAL_WALL_CELLS` (line 1030): applies boundary conditions for turbulence model
- `WALL_LOOP` (line 1175): `IW=1 to N_EXTERNAL+N_INTERNAL` — sets viscosity/thermal conductivity at wall cells
- `WALL_LOOP_SR` (line 1310): `IW=1 to N_EXTERNAL+N_INTERNAL` — Smagorinsky wall damping

### 5. VELOCITY_FLUX_KERNEL (velo_kernels.f90:325) — 11.1% of runtime

Used by: DivSetupKernelTask (predictor and corrector instances)

**Wall loops:**
- `DO IW=1,N_EXTERNAL_WALL_CELLS` (line 698): applies open boundary flux corrections

### 6. Other kernels with wall loops

- **MASS_FINITE_DIFFERENCES_NEW_KERNEL** (mass_kernels.f90:23): `WALL_LOOP_2` + `WALL_LOOP_3` for face value corrections at boundaries
- **DENSITY_KERNEL** (mass_kernels.f90:350): Wall loops + zone loops for pressure updates
- **DIVERGENCE_PART_2_KERNEL** (divg_kernels.f90:1396): Zone loops + BC_LOOP wall corrections
- **CONDENSATION_EVAPORATION_KERNEL** (fire_kernels.f90:1138): Wall loop interleaved with cell loop in species iteration
- **NO_FLUX_KERNEL** (pres_kernels.f90:737): Pure wall loop for pressure BCs

## Possible Solutions

### Approach A: Wall Cell K-Indexing

Pre-sort wall cells by their gas-cell K index. Build per-block wall cell lists at decomposition time:

```fortran
! During decomposition, build wall_cells_for_block(block_id) list
DO IW=1,N_WALL_CELLS
   KKG = BOUNDARY_COORD(WALL(IW)%BC_INDEX)%KKG
   block_id = k_to_block(KKG)
   ! Add IW to wall_cells_for_block(block_id)
ENDDO
```

Each block then processes only its own wall cells. This works IF:
- Wall cells only modify their own gas cell `(IIG,JJG,KKG)` — no cross-cell writes
- No wall cell at K=5 writes to arrays at K=35 (thin walls spanning K blocks would be a problem)
- The wall loop ordering doesn't matter (results are independent of IW iteration order)

**Applicable to:** DIVERGENCE_PART_1 wall loops (each wall cell writes to its own `(IIG,JJG,KKG)`), COMPUTE_VISCOSITY wall loops, VELOCITY_FLUX wall loops.

**NOT applicable to:** WALL_BC (processes both sides of thin walls), CONDENSATION (wall loop depends on cell loop results within species iteration), thin wall heat transfer (crosses meshes).

### Approach B: Two-Phase Split (Cell → Wall)

Split mixed kernels into:
1. **Cell phase** (block-decomposable): Pure I,J,K loops with K-range restriction
2. **Wall phase** (mesh-level): Wall loops that finalize boundary contributions

This requires restructuring the Fortran kernels so that cell loops run first, then wall loops run after block reassembly. The challenge is kernels where wall loops are interleaved with cell loops (DIVERGENCE_PART_1 has 9 wall loops mixed with cell loops).

### Approach C: Ghost-Cell Wall Overlap

Similar to MPI ghost cells: each K-block includes a small overlap zone (e.g., ±1 cell). Wall cells whose gas cell falls in the overlap zone are processed by both adjacent blocks, with the "owning" block's result taking precedence.

This adds complexity but avoids pre-sorting wall cells. However, it doesn't solve the thin-wall problem (a wall cell modifying a gas cell in a completely different K range).

### Approach D: Wall Cell Parallelism (Orthogonal to K-Blocks)

Instead of decomposing the WALL loops by K, decompose them by wall cell index:
- Split wall cells into N groups
- Each thread processes a group of wall cells

This is simpler than K-block decomposition for wall-dominant kernels but requires verifying that no two wall cells in different groups write to the same `(I,J,K)` cell. For most wall loops this holds (each wall cell writes to its unique gas cell), but thin walls and interior obstructions may create conflicts.

### Priority Targets

Based on runtime percentage, the most impactful kernels to parallelize are:

1. **DIVERGENCE_PART_1_KERNEL** (21.3%) — Most wall loops write to unique gas cells; Approach A or B could work for most of them, but the interleaved structure makes Approach B difficult.

2. **WALL_BC_PROCESS_CELLS_KERNEL** (16.2%) — Already parallelized across meshes; Approach D is most natural (each wall cell is independent except thin walls).

3. **VELOCITY_BC_PROCESS_EDGES_KERNEL** (14.5%) — Edge-based parallelism; may benefit from Approach D.

4. **COMPUTE_VISCOSITY_KERNEL + MASS_FINITE_DIFFERENCES** (11.2%) — Multiple wall loops but each writes to its own gas cell; Approach A feasible.

5. **VELOCITY_FLUX_KERNEL** (11.1%) — Single wall loop; Approach A straightforward.

Parallelizing these top 5 would cover **~74% of runtime**, potentially yielding significant speedup if the wall-cell-to-K mapping or wall-cell-index parallelism can be made safe.
