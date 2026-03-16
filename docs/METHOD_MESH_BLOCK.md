# Mesh Block Decomposition Methodology

## Goal

Add intra-mesh parallelism by decomposing individual meshes into smaller blocks that can be processed concurrently. This extends the existing inter-mesh parallelism (where each mesh is processed by one thread) to allow multiple threads to work on different regions of the same mesh.

## Context

The current Hedgehog graph processes meshes as indivisible units. Each `MeshData` token represents one full mesh, and kernel tasks process the entire mesh in one call. When the number of meshes is small relative to the number of available threads (e.g., 2 meshes on a 16-core machine), most threads sit idle.

Mesh block decomposition solves this by splitting the cell loops inside kernels into sub-ranges. Instead of `DO K=1,KBAR; DO J=1,JBAR; DO I=1,IBAR`, a block kernel processes `DO K=K1,K2; DO J=J1,J2; DO I=I1,I2`.

## Kernel Classification

Every kernel task in the graph uses one or more Fortran kernels. Each Fortran kernel falls into one of two categories:

### Mesh Block Kernel

A kernel that performs the same operation on every cell (or face) in the mesh, with no cross-cell dependencies beyond immediate neighbors. These kernels can be decomposed into blocks.

**Characteristics:**
- Loops over `(I,J,K)` with bounds `(1:IBAR, 1:JBAR, 1:KBAR)` or staggered variants `(0:IBAR, 0:JBAR, 0:KBAR)`
- Each cell's output depends only on its own value and immediate neighbors (stencil pattern)
- No global reductions across the mesh (or reductions can be done per-block and merged)
- No wall/boundary/particle structure loops interleaved with the cell loops

**Examples:** `VELOCITY_PREDICTOR_KERNEL`, `VELOCITY_CORRECTOR_KERNEL`, `PRESSURE_SOLVER_COMPUTE_RHS`

### Mesh Kernel

A kernel that must operate on the entire mesh because it loops over non-cell structures (walls, particles, zones), writes to mesh-global data, or has cross-cell dependencies that span the entire mesh.

**Characteristics:**
- Loops over walls (`DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS`)
- Loops over particles (`DO IP=1,NLP`)
- Loops over pressure zones (`DO IPZ=0,N_ZONE`)
- Writes to mesh-global scalars or boundary arrays
- Contains calls to routines that access the full mesh (e.g., FFT solvers, tridiagonal solvers)

**Examples:** `WALL_BC_PROCESS_CELLS_KERNEL`, `COMBUSTION_MODEL_KERNEL`, `PRESSURE_SOLVER_FFT`

## Step-by-Step Procedure

### Step 1: Classify the Kernels in a Task

For each kernel task in the graph, read the Fortran source of every kernel it calls. Determine whether each kernel is a mesh block kernel or a mesh kernel based on the criteria above.

Record the classification in the progress file. A task that calls multiple kernels may have a mix: some block-decomposable, some not. The task's overall classification is determined by the most restrictive kernel it contains.

**Decision matrix:**

| All kernels are block kernels | Task is a **mesh block task** — can be fully block-decomposed |
|-------------------------------|--------------------------------------------------------------|
| Some kernels are block, some are mesh | Task is a **mesh kernel task** — kept at mesh granularity unless split |
| All kernels are mesh kernels | Task is a **mesh kernel task** — no block decomposition possible |

### Step 2: Define the MeshBlockData Token

Create a new data token that represents a sub-region of a mesh:

```cpp
// Source/hedgehog/data/mesh_block_data.h
struct MeshBlockData {
    int nm;              ///< Mesh index (1-based, Fortran convention)
    int i1, i2;          ///< I-range [i1, i2] (1-based inclusive)
    int j1, j2;          ///< J-range [j1, j2] (1-based inclusive)
    int k1, k2;          ///< K-range [k1, k2] (1-based inclusive)
    double t;            ///< Current simulation time
    double dt;           ///< Current time step
    int phase;           ///< 0 = predictor, 1 = corrector

    std::shared_ptr<MeshData> originalMeshData;  ///< Parent token for reassembly
};
```

The block bounds define the sub-range of the mesh that this token represents. The kernel will process only cells within `[i1:i2, j1:j2, k1:k2]`.

**Block bounds conventions:**
- Cell-centered arrays: `i1:i2, j1:j2, k1:k2` within `1:IBAR, 1:JBAR, 1:KBAR`
- Face-centered arrays (U-faces): `i1-1:i2, j1:j2, k1:k2` (extends one cell left for staggered grid)
- Each kernel must handle the correct offset for its staggered variables

### Step 3: Create the Block Decomposition State

The decomposition state replaces the existing collector that feeds a kernel task. Instead of collecting N `MeshData` and emitting 1 `BarrierData`, it emits N*B `MeshBlockData` tokens (N meshes, B blocks per mesh):

```cpp
// Source/hedgehog/state/mesh_block_decompose_state.h
class MeshBlockDecomposeState
    : public hh::AbstractState<1, MeshData, MeshBlockData> {
public:
    MeshBlockDecomposeState(int blockSize)
        : blockSize_(blockSize) {}

    void execute(std::shared_ptr<MeshData> data) override {
        int ibar = fds_get_ibar(data->nm);
        int jbar = fds_get_jbar(data->nm);
        int kbar = fds_get_kbar(data->nm);

        // Decompose along K dimension (largest stride, best cache behavior)
        for (int k1 = 1; k1 <= kbar; k1 += blockSize_) {
            int k2 = std::min(k1 + blockSize_ - 1, kbar);
            auto block = std::make_shared<MeshBlockData>();
            block->nm = data->nm;
            block->i1 = 1; block->i2 = ibar;
            block->j1 = 1; block->j2 = jbar;
            block->k1 = k1; block->k2 = k2;
            block->t = data->t;
            block->dt = data->dt;
            block->phase = data->phase;
            block->originalMeshData = data;
            this->addResult(block);
        }
    }

private:
    int blockSize_;
};
```

### Step 4: Create the Block Reassembly State

The reassembly state collects all `MeshBlockData` tokens for a mesh and emits the original `MeshData` when all blocks for that mesh are complete:

```cpp
// Source/hedgehog/state/mesh_block_reassemble_state.h
class MeshBlockReassembleState
    : public hh::AbstractState<1, MeshBlockData, MeshData> {
public:
    MeshBlockReassembleState(int nmeshes, int blocksPerMesh)
        : nmeshes_(nmeshes), blocksPerMesh_(blocksPerMesh) {
        counts_.resize(nmeshes, 0);
        meshData_.resize(nmeshes, nullptr);
    }

    void execute(std::shared_ptr<MeshBlockData> block) override {
        int idx = block->nm - nmOffset_;
        meshData_[idx] = block->originalMeshData;
        counts_[idx]++;
        if (counts_[idx] == blocksPerMesh_) {
            counts_[idx] = 0;
            this->addResult(meshData_[idx]);
            meshData_[idx] = nullptr;
        }
    }

private:
    int nmeshes_;
    int blocksPerMesh_;
    int nmOffset_ = fds_get_lower_mesh_index();
    std::vector<int> counts_;
    std::vector<std::shared_ptr<MeshData>> meshData_;
};
```

**Note:** The number of blocks per mesh may vary if meshes have different sizes. A more robust implementation would compute blocks-per-mesh dynamically based on mesh dimensions and store the expected count per mesh.

### Step 5: Modify the Fortran Kernel

Add block-range parameters to the kernel signature. The existing full-mesh version can call the block version with full-range bounds:

```fortran
! Original kernel (full mesh)
RECURSIVE SUBROUTINE VELOCITY_PREDICTOR_KERNEL(M, DT)
TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
REAL(EB), INTENT(IN) :: DT
CALL VELOCITY_PREDICTOR_BLOCK_KERNEL(M, DT, 1, M%IBAR, 1, M%JBAR, 1, M%KBAR)
END SUBROUTINE

! Block kernel (sub-range)
RECURSIVE SUBROUTINE VELOCITY_PREDICTOR_BLOCK_KERNEL(M, DT, I1, I2, J1, J2, K1, K2)
TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
REAL(EB), INTENT(IN) :: DT
INTEGER, INTENT(IN) :: I1, I2, J1, J2, K1, K2
INTEGER :: I, J, K

DO K=K1,K2
   DO J=J1,J2
      DO I=I1-1,I2   ! Note: U-faces extend one cell left
         M%US(I,J,K) = M%U(I,J,K) - DT*( M%FVX(I,J,K) + M%RDXN(I)*(M%H(I+1,J,K)-M%H(I,J,K)) )
      ENDDO
   ENDDO
ENDDO
! ... similar for VS, WS ...
END SUBROUTINE
```

**Staggered grid handling:** FDS uses a staggered (MAC) grid where U-velocity lives on I-faces, V on J-faces, W on K-faces. When decomposing by K:
- Cell-centered variables (RHO, TMP, ZZ, H): process `K1:K2`
- W-faces: process `K1-1:K2` (one extra face at bottom of block)
- U-faces, V-faces: process `K1:K2` (same as cells)

The block kernel must handle these offsets correctly for each variable type.

### Step 6: Add C Bindings

```fortran
! fds_c_interface.f90
RECURSIVE SUBROUTINE C_FDS_VELOCITY_PREDICTOR_BLOCK_KERNEL(NM, DT, I1, I2, J1, J2, K1, K2) &
    BIND(C, NAME="fds_velocity_predictor_block_kernel")
    USE VELO_KERNELS, ONLY: VELOCITY_PREDICTOR_BLOCK_KERNEL
    USE MESH_VARIABLES, ONLY: MESHES
    INTEGER(C_INT), VALUE :: NM, I1, I2, J1, J2, K1, K2
    REAL(C_DOUBLE), VALUE :: DT
    CALL VELOCITY_PREDICTOR_BLOCK_KERNEL(MESHES(NM), DT, I1, I2, J1, J2, K1, K2)
END SUBROUTINE
```

```cpp
// fds_fortran_interface.h
void fds_velocity_predictor_block_kernel(int nm, double dt,
    int i1, int i2, int j1, int j2, int k1, int k2);
```

### Step 7: Create the Block Kernel Task

```cpp
// Source/hedgehog/task/velocity_predictor_block_kernel_task.h
class VelocityPredictorBlockKernelTask
    : public hh::AbstractTask<1, MeshBlockData, MeshBlockData> {
public:
    explicit VelocityPredictorBlockKernelTask(size_t numThreads)
        : hh::AbstractTask<1, MeshBlockData, MeshBlockData>(
              "VelPredBlockKernel", numThreads) {}

    void execute(std::shared_ptr<MeshBlockData> block) override {
        fds_velocity_predictor_block_kernel(
            block->nm, block->dt,
            block->i1, block->i2, block->j1, block->j2, block->k1, block->k2);
        this->addResult(block);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshBlockData, MeshBlockData>>
    copy() override {
        return std::make_shared<VelocityPredictorBlockKernelTask>(
            this->numberThreads());
    }
};
```

### Step 8: Wire into the Graph

Replace the existing mesh-level kernel task with the block decompose → block kernel → block reassemble pipeline:

```cpp
// Before (mesh-level):
//   prevNode -> velPredKernelTask -> nextNode

// After (block-level):
auto decompSM = std::make_shared<hh::StateManager<1, MeshData, MeshBlockData>>(
    std::make_shared<MeshBlockDecomposeState>(blockSize), "VelPredDecompose");
auto blockKernel = std::make_shared<VelocityPredictorBlockKernelTask>(kernelThreads);
auto reassembleSM = std::make_shared<hh::StateManager<1, MeshBlockData, MeshData>>(
    std::make_shared<MeshBlockReassembleState>(nmeshes, blocksPerMesh), "VelPredReassemble");

subgraph->edges(prevNode, decompSM);
subgraph->edges(decompSM, blockKernel);
subgraph->edges(blockKernel, reassembleSM);
subgraph->edges(reassembleSM, nextNode);
```

### Step 9: Build and Test

```bash
cd build_hh
cmake --build . --target fds_hh -j$(nproc)

# Run custom test suite
cd ../test_cases && python3 run_tests.py -v

# Run verification suite
python3 run_verification.py test --no-redundant --max-gold-time 30 --timeout 120 --tolerance 1e-6
```

All results must match the mesh-level baseline (no regressions).

## Decomposition Strategy

### Dimension Selection

Decompose along K (the outermost loop dimension in Fortran column-major order). This gives the best cache behavior because contiguous memory access patterns are preserved within each block.

For very flat meshes (KBAR < blockSize but IBAR or JBAR is large), consider decomposing along the largest dimension instead.

### Block Size Selection

The block size controls the granularity of parallelism:
- Too large: not enough blocks to fill all threads
- Too small: overhead from token creation and state management dominates
- Good default: `blockSize = max(1, KBAR / kernelThreads)` — creates roughly one block per thread per mesh

### Interaction with Barriers

Block decomposition happens between barriers. A barrier task (MeshExchange, HVAC, etc.) requires all meshes to be fully processed. The reassembly state before a barrier must collect all blocks for all meshes before the barrier can proceed.

The existing `CollectorState` collects N `MeshData` tokens. When block decomposition is active, the reassembly state produces N `MeshData` tokens (one per mesh, after all blocks complete), which then flow into the existing `CollectorState` unchanged.

```
MeshData(1) ─┐
MeshData(2) ─┤
             ├─> Decompose ─> BlockKernel ─> Reassemble ─> CollectorState ─> BarrierTask
MeshData(N) ─┘
```

### Handling Reductions

Some mesh block kernels compute reductions (e.g., max CFL number, sum of divergence). For block-decomposed execution:

1. Each block computes a partial reduction over its sub-range
2. The reassembly state merges partial results when collecting blocks for a mesh
3. The merged result is stored in the mesh before `MeshData` is emitted

This may require extending `MeshBlockData` with reduction fields or using per-mesh accumulators in `MESH_TYPE`.

## Relationship to Existing Patterns

This methodology builds on top of the existing parallelization:

| Level | Unit | Token | Parallelism |
|-------|------|-------|-------------|
| **Inter-mesh** (existing) | Entire mesh | `MeshData` | N meshes processed concurrently |
| **Intra-mesh** (new) | Mesh block | `MeshBlockData` | N*B blocks processed concurrently |
| **Barrier** (existing) | All meshes | `BarrierData` | Sequential cross-mesh operations |

The graph structure remains the same — barriers still collect all meshes. Only the kernel stages between barriers gain additional parallelism through block decomposition.

## Files per Conversion

Shared infrastructure (created once):
```
Source/hedgehog/data/mesh_block_data.h              Block data token
Source/hedgehog/state/mesh_block_decompose_state.h   Decompose MeshData -> MeshBlockData
Source/hedgehog/state/mesh_block_reassemble_state.h  Reassemble MeshBlockData -> MeshData
```

Per kernel conversion:
```
Source/*_kernels.f90                    Add block-range parameters to kernel
Source/hedgehog/fds_c_interface.f90     Add C binding for block kernel
Source/hedgehog/fds_fortran_interface.h Add C declaration for block kernel
Source/hedgehog/task/*_block_kernel_task.h  New block kernel task (or modify existing)
Source/hedgehog/graph/*_subgraph.h      Wire decompose -> block kernel -> reassemble
```
