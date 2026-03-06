#!/bin/bash
#
# Quick test runner for FDS Hedgehog
# Tests comparison functionality only (assumes gold files exist)
#

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMPARE_PY="$SCRIPT_DIR/compare_csv.py"
GOLD_DIR="$SCRIPT_DIR/gold"

# Test configuration
declare -A TESTS
TESTS[dancing_eddies_1mesh]="dancing_eddies_1mesh_short"
TESTS[dancing_eddies_4mesh]="dancing_eddies_4mesh_short"

# Counters
TOTAL=0
PASSED=0
FAILED=0

echo "========================================"
echo "FDS Hedgehog Quick Comparison Test"
echo "========================================"
echo

# Function to compare files
compare_files() {
    local test_name=$1
    local chid=$2
    local suffix=$3
    local gold_file="$GOLD_DIR/$test_name/${chid}${suffix}"
    local test_file="$GOLD_DIR/$test_name/${chid}${suffix}"  # Using same as gold for now

    if [ ! -f "$gold_file" ]; then
        echo -e "${YELLOW}⚠${NC} Gold file missing: ${chid}${suffix}"
        return 1
    fi

    echo -n "  Checking ${chid}${suffix}... "

    if python3 "$COMPARE_PY" "$gold_file" "$test_file" > /dev/null 2>&1; then
        echo -e "${GREEN}✓ PASS${NC}"
        return 0
    else
        echo -e "${RED}✗ FAIL${NC}"
        return 1
    fi
}

# Run tests
for test_name in "${!TESTS[@]}"; do
    chid="${TESTS[$test_name]}"

    echo "Test: $test_name"
    echo "CHID: $chid"

    test_passed=true

    # Check devc.csv
    ((TOTAL++))
    if compare_files "$test_name" "$chid" "_devc.csv"; then
        ((PASSED++))
    else
        ((FAILED++))
        test_passed=false
    fi

    # Check hrr.csv
    ((TOTAL++))
    if compare_files "$test_name" "$chid" "_hrr.csv"; then
        ((PASSED++))
    else
        ((FAILED++))
        test_passed=false
    fi

    if $test_passed; then
        echo -e "${GREEN}✓${NC} $test_name: ALL CHECKS PASSED"
    else
        echo -e "${RED}✗${NC} $test_name: SOME CHECKS FAILED"
    fi

    echo
done

# Summary
echo "========================================"
echo "Summary"
echo "========================================"
echo "Total comparisons: $TOTAL"
echo -e "Passed: ${GREEN}$PASSED${NC}"
echo -e "Failed: ${RED}$FAILED${NC}"
echo "========================================"

if [ $FAILED -eq 0 ]; then
    echo -e "${GREEN}All tests passed!${NC}"
    exit 0
else
    echo -e "${RED}Some tests failed.${NC}"
    exit 1
fi
