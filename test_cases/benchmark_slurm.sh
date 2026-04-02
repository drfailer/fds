#!/bin/bash
#SBATCH --job-name=fds-benchmark
#SBATCH --output=benchmark_%j.out
#SBATCH --error=benchmark_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --exclusive

# ── FDS Benchmark SLURM Launcher ──────────────────────────────────────────
#
# Submit on cluster:
#   sbatch benchmark_slurm.sh
#
# Or with custom options:
#   sbatch --cpus-per-task=20 benchmark_slurm.sh --omp-threads 4 8 16
#
# For multi-node MPI runs:
#   sbatch --nodes=2 --ntasks-per-node=1 benchmark_slurm.sh --mpi-procs 2
#
# Prerequisites:
#   1. Build fds_master (release) and fds_hh (release)
#   2. Run gold generation:
#      python3 test_cases/run_benchmark.py generate-gold
#   3. Make sure mpiexec is available (load MPI module if needed)
# ──────────────────────────────────────────────────────────────────────────

set -euo pipefail

# Determine script and repo directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# ── Environment Setup ────────────────────────────────────────────────────
# Uncomment/modify as needed for your cluster:

# module load gcc openmpi
# module load intel/oneapi mpi/impi

# Ensure fds_master is on PATH
if ! command -v fds_master &>/dev/null; then
    # Try common locations
    if [ -x "$REPO_ROOT/../fds-master/build/fds" ]; then
        export PATH="$REPO_ROOT/../fds-master/build:$PATH"
        # Symlink for convenience
        if [ ! -L "$REPO_ROOT/../fds-master/build/fds_master" ]; then
            ln -sf "$REPO_ROOT/../fds-master/build/fds" \
                   "$REPO_ROOT/../fds-master/build/fds_master" 2>/dev/null || true
        fi
    fi
fi

# Check executables
FDS_HH="$REPO_ROOT/build_hh/Source/hedgehog/fds_hh"
if [ ! -x "$FDS_HH" ]; then
    echo "ERROR: fds_hh not found at $FDS_HH"
    echo "Build with: cd $REPO_ROOT/build_hh && cmake --build . --target fds_hh -j\$(nproc)"
    exit 1
fi

# ── Print Environment Info ───────────────────────────────────────────────
echo "============================================================"
echo "FDS Benchmark — $(date)"
echo "============================================================"
echo "  Host:      $(hostname)"
echo "  SLURM ID:  ${SLURM_JOB_ID:-local}"
echo "  Nodes:     ${SLURM_JOB_NUM_NODES:-1}"
echo "  CPUs/task:  ${SLURM_CPUS_PER_TASK:-$(nproc)}"
echo "  fds_hh:    $FDS_HH"
echo "  fds_master: $(which fds_master 2>/dev/null || echo 'NOT FOUND')"
echo "  Python:    $(python3 --version 2>&1)"
echo "============================================================"
echo ""

# ── Detect OMP Thread Counts ────────────────────────────────────────────
# fds_master always uses MPI (nproc=nmeshes), so OMP threads are per-MPI-process.
# Avoid oversubscription: nmeshes * omp_threads should not exceed NCPUS.
# Max meshes in default benchmark suite is 8, so OMP 2 and 4 are safe up to 32 cores.
NCPUS=${SLURM_CPUS_PER_TASK:-$(nproc)}

DEFAULT_OMP_THREADS=""
for n in 2 4; do
    if [ "$n" -le "$NCPUS" ]; then
        DEFAULT_OMP_THREADS="$DEFAULT_OMP_THREADS $n"
    fi
done
DEFAULT_OMP_THREADS="${DEFAULT_OMP_THREADS# }"  # trim leading space

# ── Run Benchmark ────────────────────────────────────────────────────────
cd "$SCRIPT_DIR"

# Parse any extra arguments passed to this script (after sbatch options)
EXTRA_ARGS="$@"

# If no --omp-threads given in EXTRA_ARGS, use defaults
if [[ ! "$EXTRA_ARGS" =~ "--omp-threads" ]]; then
    EXTRA_ARGS="--omp-threads $DEFAULT_OMP_THREADS $EXTRA_ARGS"
fi

echo "Running: python3 run_benchmark.py run $EXTRA_ARGS"
echo ""

python3 run_benchmark.py run $EXTRA_ARGS

echo ""
echo "============================================================"
echo "Benchmark complete — $(date)"
echo "============================================================"

# Print report
echo ""
python3 run_benchmark.py report
