#!/bin/bash

# Test script for comparing original FDS and Hedgehog versions
# Usage: ./run_tests.sh

# Color codes for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Set paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FDS_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
FDS_ORIG="${FDS_ROOT}/build_orig/fds"
FDS_HH="${FDS_ROOT}/build_hh/fds"
TEST_DIR="${FDS_ROOT}/test_cases"
SAVED_DIR="${TEST_DIR}/saved_results"
CHID="dancing_eddies_1mesh_short"

# Create output directories
mkdir -p "${TEST_DIR}/orig_1mesh"
mkdir -p "${TEST_DIR}/orig_4mesh"
mkdir -p "${TEST_DIR}/hh_1mesh"
mkdir -p "${TEST_DIR}/hh_4mesh"

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}FDS Testing: Original vs Hedgehog${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Test 1: Original FDS - 1 mesh
echo -e "${GREEN}[1/3] Running Original FDS - 1 mesh...${NC}"
cd "${TEST_DIR}/orig_1mesh"
cp "${TEST_DIR}/${CHID}.fds" .
time mpiexec -n 1 "${FDS_ORIG}" "${CHID}.fds"
if [ $? -eq 0 ]; then
    echo -e "${GREEN}  Original FDS - 1 mesh completed${NC}"
else
    echo -e "${RED}  Original FDS - 1 mesh failed${NC}"
fi
echo ""

# Test 2: Original FDS - 4 meshes
echo -e "${GREEN}[2/3] Running Original FDS - 4 meshes...${NC}"
cd "${TEST_DIR}/orig_4mesh"
cp "${TEST_DIR}/dancing_eddies_4mesh_short.fds" .
time mpiexec -n 4 "${FDS_ORIG}" dancing_eddies_4mesh_short.fds
if [ $? -eq 0 ]; then
    echo -e "${GREEN}  Original FDS - 4 meshes completed${NC}"
else
    echo -e "${RED}  Original FDS - 4 meshes failed${NC}"
fi
echo ""

# Test 3: Hedgehog FDS - 1 mesh
echo -e "${GREEN}[3/3] Running Hedgehog FDS - 1 mesh...${NC}"
cd "${TEST_DIR}/hh_1mesh"
cp "${TEST_DIR}/${CHID}.fds" .
time mpiexec -n 1 "${FDS_HH}" "${CHID}.fds"
if [ $? -eq 0 ]; then
    echo -e "${GREEN}  Hedgehog FDS - 1 mesh completed${NC}"
else
    echo -e "${RED}  Hedgehog FDS - 1 mesh failed${NC}"
fi
echo ""

# Compare results against saved baselines
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Comparing Against Saved Results${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

PASS=0
FAIL=0

compare_files() {
    local label="$1"
    local test_dir="$2"
    local saved_dir="$3"

    if [ ! -d "$saved_dir" ]; then
        echo -e "${RED}  No saved results in ${saved_dir}${NC}"
        return
    fi

    # Compare CSV files (devc, hrr) — these should be byte-identical
    for suffix in devc.csv hrr.csv; do
        local test_file="${test_dir}/${CHID}_${suffix}"
        local saved_file="${saved_dir}/${CHID}_${suffix}"
        if [ ! -f "$test_file" ]; then
            continue
        fi
        if diff -q "$test_file" "$saved_file" > /dev/null 2>&1; then
            echo -e "${GREEN}  ${label} ${suffix}: IDENTICAL${NC}"
            PASS=$((PASS + 1))
        else
            echo -e "${RED}  ${label} ${suffix}: DIFFERS${NC}"
            FAIL=$((FAIL + 1))
        fi
    done

    # Compare binary slice files
    for saved_file in "${saved_dir}"/${CHID}_*.sf; do
        [ -f "$saved_file" ] || continue
        local base
        base=$(basename "$saved_file")
        local test_file="${test_dir}/${base}"
        if [ ! -f "$test_file" ]; then
            continue
        fi
        if diff -q "$test_file" "$saved_file" > /dev/null 2>&1; then
            echo -e "${GREEN}  ${label} ${base}: IDENTICAL${NC}"
            PASS=$((PASS + 1))
        else
            echo -e "${RED}  ${label} ${base}: DIFFERS${NC}"
            FAIL=$((FAIL + 1))
        fi
    done
}

compare_files "orig_1mesh" "${TEST_DIR}/orig_1mesh" "${SAVED_DIR}/orig_1mesh"
compare_files "hh_1mesh"   "${TEST_DIR}/hh_1mesh"   "${SAVED_DIR}/hh_1mesh"

echo ""
echo -e "${BLUE}========================================${NC}"
if [ "$FAIL" -eq 0 ]; then
    echo -e "${GREEN}All ${PASS} comparisons passed.${NC}"
else
    echo -e "${RED}${FAIL} comparison(s) failed, ${PASS} passed.${NC}"
fi
echo -e "${BLUE}========================================${NC}"
