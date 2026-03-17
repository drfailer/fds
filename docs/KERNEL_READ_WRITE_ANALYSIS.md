# Kernel Read/Write Analysis for K-Block Decomposition

## Key Insight

The original WALL_LOOP_BARRIER.md classified most kernels as "Mesh" because they
contain wall loops. However, a deeper analysis of **what** each wall loop reads and
writes reveals that many wall loops are actually K-decomposable:

1. **Ghost cell writes** — Write to exterior boundary cells (K=0, K=KBP1). Safe for
   K-decomposition: ghost cells are outside interior K ranges.
2. **Face-value writes** — Write to face arrays (FX, FY, KDTDX, etc.) at the face
   adjacent to a wall's gas cell. Location is deterministic: (IIG-1, JJG, KKG) for
   IOR=1, etc. Can be K-partitioned by gas cell K.
3. **Gas cell writes** — Write to field arrays at (IIG, JJG, KKG). Can be
   K-partitioned by KKG.
4. **Gas cell accumulations** — `+=` to field arrays at (IIG, JJG, KKG). K-safe if
   each gas cell's wall contributions all fall in the same K block (true unless thin
   walls span K blocks — rare).
5. **Wall cell property writes** — Write to B1%property (per-wall-cell, not per-grid).
   Always safe for any decomposition.
6. **Off-wall corrections** — Write to face arrays at (II+2, KK-2, etc.). Can cross
   K-block boundaries. Needs ghost zones or sequential post-processing.
7. **Global/zone writes** — CONNECTED_ZONES, D_SUM_LOC, PBAR, etc. Must stay
   sequential or use per-block accumulators.

Cross-cell **reads** (e.g., `UU(I+1,J,K)`, `WW(I,J,K-1)`) are never a problem for
K-decomposition — they only require read access to neighboring cells, which is
guaranteed by the staggered grid layout.

---

## 1. VELOCITY_FLUX_KERNEL (velo_kernels.f90:325-913) — 11.1% runtime

**Used by:** DivSetupKernelTask (predictor + corrector)

### Cell loops (main body):

| Loop | Range | Reads | Writes | Type |
|------|-------|-------|--------|------|
| Vorticity/stress (367) | K=0:KBAR | UU,VV,WW at neighbors; MU at neighbors | OMX,OMY,OMZ,TXY,TXZ,TYZ (WORK1-6) at (I,J,K) | In-cell |
| FVX (406) | K=1:KBAR | OMY,OMZ,TXY,TXZ at K,K-1; EDGE%OMEGA/TAU; RHOP,DP,VV,WW at neighbors | FVX(I,J,K) | In-cell |
| FVY (462) | K=1:KBAR | OMX,OMZ,TYZ,TXY at K,K-1; EDGE data; RHOP,DP,UU,WW at neighbors | FVY(I,J,K) | In-cell |
| FVZ (518) | K=0:KBAR | OMX,OMY,TXZ,TYZ at neighbors; EDGE data; RHOP,DP,UU,VV at neighbors | FVZ(I,J,K) | In-cell |

### CONTAINS subroutines:

| Subroutine | Condition | Wall loops | Write pattern |
|------------|-----------|------------|---------------|
| DIRECT_FORCE (585) | `ABS(FVEC)>0` | None | FVX/FVY/FVZ at (I,J,K) — in-cell |
| CORIOLIS_FORCE (671) | `ABS(OVEC)>0` | 1 wall loop (698): `IW=1,N_EXTERNAL` | UP/VP/WP (WORK7-9) at **ghost cells** (BC%II,JJ,KK) |
| PATCH_VELOCITY_FLUX (776) | `PATCH_VELOCITY` | None | FVX/FVY/FVZ within device K1:K2 ranges |
| MMS_VELOCITY_FLUX (745) | `PERIODIC_TEST==7` | None | FVX/FVZ at (I,J,K) — in-cell |

### Assessment: **BLOCK-DECOMPOSABLE (Two-Phase Split)**

The main body is 100% cell loops with in-cell writes. The only wall loop is in
CORIOLIS_FORCE (rare: rotating reference frame only), and it writes to **ghost cells
only** — the WORK arrays at exterior boundary positions. These ghost values are then
read by subsequent cell loops.

**Split strategy:**
- Phase 1 (sequential): CORIOLIS_FORCE wall loop (if active) — set ghost cells
- Phase 2 (K-parallel): All cell loops (vorticity, FVX, FVY, FVZ, forces)

For most simulations (no Coriolis), Phase 1 is empty and the entire kernel is parallel.

The EDGE reads (`M%EDGE(IEYP)%OMEGA`) access edge data by cell-based index lookup —
the edge index is determined by CELL_INDEX(I,J,K), so each cell reads its own edges.
This is safe for K-decomposition.

---

## 2. COMPUTE_VISCOSITY_KERNEL (velo_kernels.f90:921-1273) — 11.2% runtime

**Used by:** PredStep1KernelTask, CorrStep1KernelTask

### Cell loops:

| Loop | Range | Writes | Type |
|------|-------|--------|------|
| MU_DNS (970) | K=1:KBAR | MU_DNS(I,J,K) | In-cell |
| STRAIN_RATE (1283) | K=1:KBAR | STRAIN_RATE(I,J,K) | In-cell |
| Turb model (998-1152) | K=1:KBAR | MU(I,J,K) | In-cell |
| KRES (1157) | K=1:KBAR | KRES(I,J,K) | In-cell |

### Wall loops:

| Loop | Range | Writes | Write type |
|------|-------|--------|------------|
| Deardorff ghost (1030) | `IW=1,N_EXTERNAL` | UP/VP/WP at (BC%II,JJ,KK) | **Ghost cell** |
| WALL_LOOP (1175) | `IW=1,N_EXT+N_INT` | MU(IIG,JJG,KKG), KRES(II,JJ,KK), MU(II,JJ,KK) | **Gas cell** + ghost mirror |
| WALL_LOOP_SR (1310) | `IW=1,N_EXT+N_INT` | STRAIN_RATE(IIG,JJG,KKG) | **Gas cell** |
| Corner mirroring (1248) | Array sections | MU, KRES at edges/corners | **Mesh boundary** |

### Assessment: **BLOCK-DECOMPOSABLE (Two-Phase Split)**

**Split strategy:**
- Phase 1 (K-parallel): MU_DNS, STRAIN_RATE, turb model, KRES cell loops
- Phase 2 (sequential): WALL_LOOP (overwrite MU/KRES at boundary gas cells), corner
  mirroring

WALL_LOOP writes to gas cells at (IIG,JJG,KKG) — the cell adjacent to the wall. This
is an overwrite (not accumulation), using CELL_COUNTER for weighted averaging when
multiple walls share a gas cell. Running this sequentially after the parallel phase is
correct because the wall loop replaces the cell-loop value.

**Complication:** Deardorff model (the default!) has a wall loop (line 1030) between
the UP/VP/WP cell loop and TEST_FILTER. This loop writes ghost cells only, so it could
either:
- (a) Be run as a pre-processing step before K-decomposition
- (b) Be included in each block if blocks extend to ghost cells

The FILL_EDGES_KERNEL and TEST_FILTER_KERNEL calls (1048-1061) operate on full mesh
arrays. These would need to stay in the sequential phase or be adapted.

**Practical split for non-Deardorff (CONSMAG, VREMAN, WALE):**
- Phase 1 (parallel): All cell loops
- Phase 2 (sequential): Wall loops + corner mirroring
- Clean separation, no interleaving

**For Deardorff:**
- More complex due to wall loop → FILL_EDGES → TEST_FILTER interleaving in step 3
- Would require splitting the Deardorff block into: cell loop → [barrier] → wall+filter → [barrier] → cell loop
- May not be worth the complexity

---

## 3. MASS_FINITE_DIFFERENCES_NEW_KERNEL (mass_kernels.f90:23-341) — part of Step1

### Cell loops:

| Loop | Range | Writes | Type |
|------|-------|--------|------|
| RHO_Z_P (69) | extended range | WORK_PAD(I,J,K) | In-cell |
| GET_SCALAR_FACE_VALUE (82-84) | structured | FX/FY/FZ(:,:,:,N) | Face value |
| MW correction (315-337) | K=0:KBAR | FX/FY/FZ(:,:,:,N) | In-cell |

### Wall loops:

| Loop | Range | Writes | Write type |
|------|-------|--------|------------|
| WALL_LOOP_2 (89) | `IW=1,N_EXT+N_INT` | FX/FY/FZ at (IIG±1,JJG±1,KKG±1,N) | **Face at gas cell** |
| | | FX/FY/FZ at (II±2,JJ±2,KK±2,N) — off-wall | **Off-wall face** |
| WALL_LOOP_3 (214) | `IW=1,N_EXT+N_INT` | Same pattern for FX/FY/FZ(:,:,:,0) | **Face at gas cell** + off-wall |

### Assessment: **MIXED — Split possible but off-wall corrections complicate**

The main face value writes are at the wall's gas cell face — K-partitionable by KKG.
However, the **off-wall corrections** write to positions like (II+1,JJ,KK) or
(II,JJ,KK-2) relative to the ghost cell. For K-oriented walls (IOR=±3), the off-wall
write at KK±2 could cross a K-block boundary.

**Mitigation:** The off-wall corrections are conditional (`UU(II+1,JJ,KK)>0` etc.) and
only affect cells immediately adjacent to walls. A 1-cell ghost zone in K would handle
this, or the off-wall corrections could be post-processed sequentially.

**Split strategy:**
- Phase 1 (K-parallel): RHO_Z_P computation + GET_SCALAR_FACE_VALUE
- Phase 2 (sequential or K-partitioned with ghost): Wall face corrections + off-wall corrections
- Phase 3 (K-parallel): MW correction loop

---

## 4. DENSITY_KERNEL (mass_kernels.f90:350-694) — 1.2% runtime

### Cell loops:

| Loop | Range | Writes | Type |
|------|-------|--------|------|
| Species advection (422-435) | K=1:KBAR | ZZS(I,J,K,N) or ZZ(I,J,K,N) | In-cell |
| RHOS/RHO sum (481-488) | K=1:KBAR | RHOS(I,J,K) or RHO(I,J,K) | In-cell |
| ZZ normalization (498-505) | K=1:KBAR | ZZS(I,J,K,:)/RHOS or ZZ/RHO | In-cell |
| RSUM (519-527) | K=1:KBAR | RSUM(I,J,K) | In-cell |
| TMP from EOS (531-538) | K=1:KBAR | TMP(I,J,K) | In-cell |

### Wall loops:

| Loop | Range | Writes | Write type |
|------|-------|--------|------------|
| WALL_LOOP (406) | `IW=1,N_EXTERNAL` | UU/VV/WW WORK at (IIG±1,JJG±1,KKG±1) | **Face velocity at gas cell** |

### Other:

| Operation | Writes | Type |
|-----------|--------|------|
| PBAR_S update (513-515) | PBAR_S(:,IPZ) | **Zone-level** (1D, not per-cell) |
| CHECK_MASS_DENSITY | CLIP_RHOMIN/RHOMAX mesh flags | **Mesh-level reduction** |

### Assessment: **BLOCK-DECOMPOSABLE (Two-Phase Split)**

All cell loops are in-cell writes. The wall loop only affects velocity WORK arrays at
boundary faces (INTERPOLATED_BOUNDARY cells only). CHECK_MASS_DENSITY is a mesh-level
reduction (min/max scan over all cells).

**Split strategy:**
- Phase 1 (sequential): Wall loop (set boundary velocities), PBAR update
- Phase 2 (K-parallel): All cell loops (advection, sum, normalize, RSUM, TMP)
- Phase 3 (sequential): CHECK_MASS_DENSITY (mesh-level flags)

---

## 5. DIVERGENCE_PART_1_KERNEL (divg_kernels.f90:27-1387) — 21.3% runtime

This is the most complex kernel with 9+ wall loops interleaved with cell loops.

### Wall loops — detailed write analysis:

| # | Loop | Range | Key writes | Write type | K-safe? |
|---|------|-------|------------|------------|---------|
| 1 | WALL_LOOP3 (99) | `IW=1,N_EXT+N_INT` | B1%U_NORMAL_S, B1%U_NORMAL | Wall cell property | Yes |
| 2 | WALL_LOOP (195) | `IW=1,N_EXT+N_INT` | RHO_D_DZDX at (IIG-1,JJG,KKG), etc. | Face at gas cell | Yes |
| 3 | WALL_LOOP_2 (315) | `IW=1,N_EXT+N_INT` | RHO_D_DZDX, H_RHO_D_DZDX at faces; B1%RHO_D_DZDN_F | Face + wall prop | Yes |
| 4 | BOUNDARY_LOOP (469) | `IW=1,N_EXTERNAL` | KP at (BC%II,JJ,KK) | Ghost cell | Yes |
| 5 | CORRECTION_LOOP (515) | `IW=1,N_EXT+N_INT` | **DP(IIG,JJG,KKG)** `+=`; KDTDX at faces | **Gas cell accum** | Mostly* |
| 6 | WALL_LOOP (807) in ENTHALPY_ADVECTION_NEW | `IW=1,N_EXT+N_INT` | FX_H_S off-wall; **U_DOT_DEL_RHO_H_S(IIG,JJG,KKG)** `+=` | Off-wall + **gas cell accum** | Mostly* |
| 7 | WALL_LOOP_2 (988) in SPECIES_ADVECTION_PART_1_NEW | `IW=1,N_EXT+N_INT` | FX_ZZ/FY_ZZ/FZ_ZZ at ghost+off-wall | Face + off-wall | Mostly* |
| 8 | WALL_LOOP (1195) in SPECIES_ADVECTION_PART_2 | `IW=1,N_EXT+N_INT` | **U_DOT_DEL_RHO_Z(IIG,JJG,KKG)** `+=` | **Gas cell accum** | Mostly* |
| 9 | WALL_LOOP4 (734) | `IW=1,N_EXT+N_INT` | D_SUM_LOC(IPZ), P_SUM_LOC(IPZ), U_SUM_LOC(IPZ) | **Zone accum** | Per-block accum |

*"Mostly" = safe if wall cells are K-partitioned by their gas cell's KKG. Thin walls
spanning K blocks are the exception.

### Key problems:

1. **Gas cell accumulations** (#5, #6, #8): Multiple wall cells can `+=` to the same
   gas cell. If all wall cells at a given (IIG,JJG,KKG) have the same KKG (which they
   do — the gas cell has a fixed K), this is safe for K-partitioning.

2. **Off-wall corrections** (#6, #7): Write to face arrays at (II±2, KK±2) which can
   cross block boundaries. Same issue as MASS_FINITE_DIFFERENCES.

3. **Interleaving**: The wall loops are not at the start/end — they occur between cell
   loops within SPECIES_LOOP and SPECIES_ADVECTION iterations.

4. **Zone accumulators** (#9): Already solved with per-mesh local arrays. Could extend
   to per-block local arrays.

5. **MERGE_PRESSURE_ZONES**: Writes to global CONNECTED_ZONES. Must stay sequential.

### Assessment: **VERY DIFFICULT to split, but cell-loop portions are decomposable**

The interleaving of wall loops and cell loops within the species iteration makes a
clean split impractical. However, the individual sub-operations (diffusive flux
computation, enthalpy advection, species advection) follow a pattern of:
1. Cell loop computes face values
2. Wall loop corrects face values at boundaries
3. Cell loop uses corrected face values

Each such pair could be split, but the number of synchronization points would be high.

---

## 6. VELOCITY_BC_PROCESS_EDGES_KERNEL (velo_kernels.f90:1544) — 14.5% runtime

### Structure: Single loop over edges (IE=1, EDGE_COUNT)

| Operation | Writes | Type |
|-----------|--------|------|
| EDGE%OMEGA/TAU | Edge data structure | Per-edge |
| UU/VV/WW at boundaries | Velocity at K=0, K=KBP1, etc. | Boundary face |

### Assessment: **MESH-LEVEL (edge-indexed, not K-decomposable)**

The edge loop iterates over a flat edge index, not (I,J,K). Each edge maps to
specific (II,JJ,KK) coordinates, but the edge index doesn't correlate with K.
The writes to UU/VV/WW are at domain boundaries.

Could theoretically be K-partitioned by edge coordinate (each edge has II,JJ,KK),
but the complex branching logic and wall cell lookups make this impractical.

---

## 7. WALL_BC_PROCESS_CELLS_KERNEL — 16.2% runtime

### Structure: Single loop over wall cells (IW=1, N_WALL_CELLS)

Each wall cell is processed independently (heat transfer, species boundary conditions).
Writes are to wall cell properties (B1, B2) and to gas cell arrays at (IIG,JJG,KKG).

### Assessment: **K-DECOMPOSABLE (Approach A: Wall Cell K-Indexing)**

Each wall cell writes to its own gas cell at (IIG,JJG,KKG). If wall cells are
pre-sorted by KKG into per-block lists, each block processes its own wall cells.

**Complication:** Thin walls — a single physical wall generates two wall cells on
opposite sides, both mapping to the same gas cell. These must be in the same K block
(which they are, since both map to the same KKG).

---

## 8. DIVERGENCE_PART_2_KERNEL (divg_kernels.f90:1396) — 2.3% runtime

### Assessment: **MESH-LEVEL**

Contains zone loops (D_PBAR_DT updates over pressure zones) that must run before cell
loops. Small runtime — not worth splitting.

---

## 9. PRESSURE_SOLVER (pres_kernels.f90) — 4.8% runtime

### Assessment: **MESH-LEVEL**

FFT/ULMAT solvers operate on the full mesh. Already parallelized via pressure iteration
sub-graph.

---

## 10. CONDENSATION_EVAPORATION_KERNEL (fire_kernels.f90:1138)

### Assessment: **MESH-LEVEL**

Wall loop interleaved with cell loop within species iteration, sharing interim arrays.
Cannot cleanly separate.

---

## Summary: Revised Classification

| # | Kernel | % Runtime | Old Class | New Class | Split Strategy |
|---|--------|-----------|-----------|-----------|----------------|
| 1 | VELOCITY_FLUX | 11.1% | Mesh | **Block (2-phase)** | Ghost wall loop → parallel cell loops |
| 2 | COMPUTE_VISCOSITY | 11.2% | Mesh | **Block (2-phase)** | Parallel cell loops → sequential wall fixup |
| 3 | MASS_FINITE_DIFFS | (part of Step1) | Mesh | **Mixed (2-phase)** | Parallel face values → sequential wall corrections |
| 4 | DENSITY | 1.2% | Mesh | **Block (2-phase)** | Sequential wall+PBAR → parallel cells → sequential check |
| 5 | DIVERGENCE_PART_1 | 21.3% | Mesh | **Mesh** | Too many interleaved wall/cell loops |
| 6 | VELOCITY_BC_EDGES | 14.5% | Mesh | **Mesh** | Edge-indexed, not K-decomposable |
| 7 | WALL_BC | 16.2% | Mesh | **Block (Approach A)** | Pre-sort wall cells by KKG, each block processes its own |
| 8 | DIV_PART_2 | 2.3% | Mesh | Mesh | Zone loops, small runtime |
| 9 | PRESSURE_SOLVER | 4.8% | Mesh | Mesh | FFT/ULMAT global |
| 10 | CONDENSATION | (small) | Mesh | Mesh | Interleaved wall+cell |

### Newly decomposable runtime:

```
VELOCITY_FLUX:       11.1% (clean 2-phase split)
COMPUTE_VISCOSITY:   11.2% (2-phase, Deardorff needs care)
WALL_BC:             16.2% (wall cell K-indexing)
DENSITY:              1.2% (2-phase split)
MASS_FINITE_DIFFS:   ~3%   (2-phase, off-wall complication)
─────────────────────────
TOTAL:               ~42.7%
```

Combined with the already-decomposed velocity predictor/corrector (0.2%) and the
remaining mesh-level kernels (DIVERGENCE_PART_1 at 21.3%, VELOCITY_BC_EDGES at 14.5%),
**~43% of runtime could be K-decomposed** vs the previous 0.2%.

### Remaining mesh-level bottlenecks:

```
DIVERGENCE_PART_1:   21.3% — interleaved wall/cell loops (hardest to split)
VELOCITY_BC_EDGES:   14.5% — edge-indexed loops
PRESSURE_SOLVER:      4.8% — FFT/ULMAT (already parallel via sub-graph)
DIV_PART_2:           2.3% — zone-level operations
Other:                ~7%  — I/O, radiation, combustion, particles
```

### Priority implementation order:

1. **WALL_BC** (16.2%) — Approach A (wall cell K-indexing) is cleanest
2. **VELOCITY_FLUX** (11.1%) — Main body has zero wall loops; Coriolis is rare
3. **COMPUTE_VISCOSITY** (11.2%) — Non-Deardorff models have clean 2-phase split
4. **DENSITY** (1.2%) — Simple split, low impact
5. **MASS_FINITE_DIFFS** (~3%) — Off-wall complications, moderate impact
