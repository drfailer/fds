# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Fire Dynamics Simulator (FDS) is a large-eddy simulation (LES) code for low-speed flows, with an emphasis on smoke and heat transport from fires. The codebase is written in Fortran 2018 and uses MPI for parallelization across multiple meshes.

## Build System

FDS supports both traditional Makefile and CMake build systems.

### Makefile Build

Navigate to a build directory under `Build/` and run the build script:

```bash
cd Build/impi_intel_linux_openmp_db
./make_fds.sh
```

Build directory names indicate: MPI implementation (impi/ompi), compiler (intel/gnu), OS (linux/win/osx), and mode:
- No suffix = release (optimized)
- `_db` = debug mode
- `_dv` = development (low optimization)
- `_openmp` = compiled with OpenMP directives

Common build targets from `Build/makefile`:
- `make impi_intel_linux` - Intel MPI + Intel compiler, Linux, release
- `make impi_intel_linux_db` - Debug version
- `make ompi_gnu_linux` - OpenMPI + GNU compiler, Linux

### CMake Build

```bash
cmake --preset default
cmake --build --preset default
```

CMake options:
- `USE_HYPRE=ON/OFF` - Use hypre library (default: ON)
- `USE_SUNDIALS=ON/OFF` - Use sundials library (default: ON)
- `USE_OPENMP=ON/OFF` - Enable OpenMP (default: ON)
- `USE_SYSTEM_HYPRE=ON/OFF` - Use system hypre vs download
- `USE_SYSTEM_SUNDIALS=ON/OFF` - Use system sundials vs download

## Code Architecture

### Source Code Structure

All Fortran source is in `Source/`. The main entry point is `main.f90`. The code is organized into functional modules:

**Core Data Structures:**
- `type.f90` (TYPES) - All derived types (MESH_TYPE, WALL_TYPE, LAGRANGIAN_PARTICLE_TYPE, etc.)
- `mesh.f90` (MESH_VARIABLES) - MESH_TYPE definition containing all per-mesh variables
- `data.f90` (OUTPUT_DATA) - Output quantity definitions
- `cons.f90` (GLOBAL_CONSTANTS) - Physical constants and flags

**Numerical Solvers:**
- `pres.f90` - Pressure solvers (FFT, ULMAT, GLMAT)
- `divg.f90` - Divergence calculations
- `velo.f90` - Velocity predictor
- `mass.f90` - Density and species advection
- `pois.f90` - Poisson equation solver

**Physics Modules:**
- `fire.f90` - Combustion and heat release
- `radi.f90` - Radiation transport
- `turb.f90` - Turbulence modeling (Smagorinsky LES)
- `wall.f90` - Boundary layer and wall heat transfer
- `part.f90` - Lagrangian particle tracking
- `chem.f90` - Chemistry integration
- `hvac.f90` - HVAC network solver

**I/O and Setup:**
- `read.f90` - FDS input file parsing
- `init.f90` - Mesh initialization
- `dump.f90` - Output and restart files

**Complex Geometry:**
- `ccib.f90` - Cut-cell Cartesian immersed boundary method
- `geom.f90` - Geometry intersection calculations

### Program Flow

1. **Initialization** (main.f90:80-556):
   - MPI/OpenMP setup
   - Input file parsing (`READ_DATA`)
   - Mesh allocation and coordinate setup
   - Wall/boundary initialization
   - Device, particle, radiation initialization

2. **Main Time-Stepping Loop** (main.f90:561+):
   - **PREDICTOR step**: Estimates state at n+1
     - Particle insertion and mass transport
     - Species density prediction
     - Velocity flux calculations
     - Wall boundary conditions
     - Divergence calculations
     - **Pressure iteration**: Solves Poisson equation iteratively until velocity error < VELOCITY_TOLERANCE
     - Velocity prediction with CFL-based time step control

   - **CORRECTOR step**: Refines solution
     - Species mass corrections
     - Energy equations
     - Particle movement and interactions
     - Output/diagnostics

3. **Finalization**: Diagnostic dumps, restart file writing

### Key Architectural Patterns

**Mesh-Based Parallelization:**
- Each MPI process owns a subset of meshes (NMESHES)
- Primary data structure: `MESH_TYPE` (MESHES array indexed by NM)
- Loop pattern: `DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX`
- Boundary data exchanged via `MESH_EXCHANGE(CODE)` where CODE specifies what to exchange (1=species, 3=velocity/pressure, 5=momentum fluxes, 7=particles, etc.)
- Neighboring mesh data stored in `OMESH_TYPE` arrays

**Predictor-Corrector Time Stepping:**
- Two-stage approach for implicit time integration
- Pressure solved iteratively within predictor if velocity error exceeds tolerance
- Time step adapted based on CFL condition checked in `VELOCITY_PREDICTOR`

**Pressure Solver Architecture:**
- Three solver strategies (selected by PRES_FLAG):
  - FFT_FLAG: Fast Fourier Transform (periodic/structured)
  - ULMAT_FLAG: Unstructured Local Matrix (sparse, general domains)
  - GLMAT_FLAG: Global Matrix (across all meshes)
- Each pressure zone may have separate matrix
- MUNKH array maps cell (I,J,K) to pressure unknown index

**Wall Cell Boundary Treatment:**
- `WALL_TYPE` contains boundary data for all wall cells on a mesh
- BOUNDARY_COORD: Indices for ghost cell and gas cell, orientation
- BOUNDARY_ONE_D: 1-D solid-phase conduction/pyrolysis solver for thermally thick surfaces
- Wall arrays dynamically grow as needed (N_WALL_CELLS vs N_WALL_CELLS_DIM)

**Lagrangian Particle Tracking:**
- Particles advected through velocity field with two-way coupling
- Particle-gas interaction via drag forces (FVX_D, FVY_D, FVZ_D) and heat/mass transfer (M_DOT_PPP, Q)
- Workflow: `INSERT_ALL_PARTICLES` → `MOVE_PARTICLES` → `PARTICLE_MASS_ENERGY_TRANSFER` → `REMOVE_PARTICLES`

**Complex Geometry (CC_IBM):**
- Cut-cell immersed boundary method for arbitrary geometry on Cartesian grids
- Special treatment in Poisson solver and advection schemes
- Activated when CC_IBM=.TRUE.

## Testing and Verification

### Quick Debug Check

Before making substantial changes, run verification cases in debug mode for a few time steps:

```bash
# Compile debug versions
cd Build/impi_intel_linux_db
./make_fds.sh

# Clean verification directory
cd ../../Verification
git clean -dxf  # WARNING: Erases all uncommitted files

# Run cases for 2 time steps in debug mode
cd scripts
./Run_FDS_Cases.sh -m 2 -d -q firebot

# Check for Fortran runtime errors
cd ..
grep forrtl */*err
```

### Full Verification Suite

```bash
# Compile release versions
cd Build/impi_intel_linux_openmp
./make_fds.sh

# Clean and run full suite
cd ../../Verification
git clean -dxf
cd scripts
./Run_FDS_Cases.sh -q firebot  # Takes a few hours

# Test restart feature
./Run_FDS_Cases.sh -r

# Process results with Matlab
cd ../../Utilities/Matlab
# Run FDS_verification_script.m in Matlab (10-15 minutes)

# Check for failures
cd ../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/Scatterplots
# Open verification_scatterplot_output.csv
```

### Validation Cases

Validation cases are in `Validation/` with sub-folders for each test series. Each contains:
- `Run_All.sh` - Creates `Current_Results/` and runs FDS input files
- `Process_Output.sh` - Copies output to separate repository

Run validation cases:
```bash
cd Validation
./Run_Serial.sh      # Single-process jobs (run first)
./Run_Parallel.sh    # Multi-process jobs (run after serial)
./Process_All_Output.sh  # Process completed cases
```

## Development Guidelines

### Code Modification Patterns

**When modifying source code:**
1. All changes must successfully execute the full V&V suite
2. Avoid large speed and memory penalties
3. Work within the existing FDS framework
4. Maintain broad applicability and ease of use
5. Document new algorithms
6. Be prepared to support issue resolution for 6-12 months after release

**Likely to be accepted:**
- New DEVC or CTRL types that don't break existing features
- Minor changes leveraging existing code (e.g., making constant parameters RAMP-dependent)
- Performance improvements without changing fundamental approach
- New physical submodels that meet restrictions and demonstrate improvement

**Unlikely to be accepted:**
- Changes to fundamental structure (e.g., compressible flow solver)
- Links to specialized precompiled or non-public domain libraries

### Git Workflow

FDS uses a standard git workflow. When making changes:
1. Create a feature branch from master
2. Make focused commits
3. Test with verification suite
4. Submit pull request to master branch

## Additional Resources

- [FDS-SMV Website](https://pages.nist.gov/fds-smv/)
- [Firebot Build Status](https://pages.nist.gov/fds-smv/firebot_status.html)
- [Discussions](https://github.com/firemodels/fds/discussions)
- [FDS Issues](https://github.com/firemodels/fds/issues)
- [Developer Commit Guidelines](https://github.com/firemodels/fds/wiki/Developer-Commit-Guidelines)
- [FDS Road Map](https://github.com/firemodels/fds/wiki/FDS-Road-Map)
- [FDS Verification Process](https://github.com/firemodels/fds/wiki/FDS-Verification-Process)
- [FDS Validation Process](https://github.com/firemodels/fds/wiki/FDS-Validation-Process)
