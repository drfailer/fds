# FDS Hedgehog Test System

Automated testing framework for FDS with Hedgehog integration.

## Directory Structure

```
test_cases/
├── inputs/              # Test input files (.fds)
├── gold/                # Gold standard output files (reference)
│   ├── test_name_1/     # Gold files for test 1
│   └── test_name_2/     # Gold files for test 2
├── run/                 # Test execution directory (temporary)
├── compare_csv.py       # CSV comparison utility
└── run_tests.py         # Main test runner
```

## Test Cases

| Test Name | Meshes | Time | Description |
|-----------|--------|------|-------------|
| `dancing_eddies_1mesh` | 1 | 0.1s | Single mesh Dancing Eddies |
| `dancing_eddies_2mesh` | 2 | 2.0s | Embedded mesh test |
| `multiple_reac_3mesh` | 3 | 0.0005s | Multiple chemical reactions |
| `dancing_eddies_4mesh` | 4 | 0.1s | Four mesh Dancing Eddies |
| `species_props_5mesh` | 5 | 0.05s | Species properties |

## Usage

### Running Tests

Run all tests:
```bash
cd test_cases
./run_tests.py
```

Run specific test:
```bash
./run_tests.py --test dancing_eddies_4mesh
```

Run multiple specific tests:
```bash
./run_tests.py --test dancing_eddies_1mesh --test dancing_eddies_4mesh
```

Verbose output:
```bash
./run_tests.py --verbose
```

Custom tolerance:
```bash
./run_tests.py --tolerance 1e-12
```

### Generating Gold Files

Generate gold files using original FDS (from master branch):
```bash
# First, ensure you have built the original FDS
cd Build/ompi_gnu_linux_db
./make_fds.sh

# Then generate gold files
cd ../../test_cases
./run_tests.py --generate-gold
```

Generate gold files using fds_hh (hedgehog version):
```bash
./run_tests.py --generate-gold --use-hedgehog-gold
```

Generate gold for specific test only:
```bash
./run_tests.py --generate-gold --test dancing_eddies_1mesh
```

## CSV Comparison Tool

The `compare_csv.py` script compares FDS CSV outputs while ignoring timing information.

Usage:
```bash
./compare_csv.py file1.csv file2.csv
./compare_csv.py file1.csv file2.csv --tolerance 1e-12
./compare_csv.py file1.csv file2.csv --verbose
```

Features:
- Ignores timing columns (Time, Step, CPU, etc.)
- Supports custom numerical tolerance
- Detects both absolute and relative differences
- Reports first 10 failures for debugging

## Adding New Tests

1. **Create input file** in `inputs/` directory
2. **Add test configuration** to `run_tests.py`:
```python
TEST_CASES = {
    'my_new_test': {
        'input': 'my_test.fds',
        'meshes': 4,
        'description': 'Description of test',
        'compare_files': ['_devc.csv', '_hrr.csv']  # Files to compare
    },
}
```
3. **Generate gold files**:
```bash
./run_tests.py --generate-gold --test my_new_test
```
4. **Run test**:
```bash
./run_tests.py --test my_new_test
```

## Expected Output Comparison

The test runner compares these output files (when present):
- `*_devc.csv` - Device outputs (thermocouples, velocity measurements, etc.)
- `*_hrr.csv` - Heat release rate
- `*_steps.csv` - Time step information
- Additional CSV files can be added per test

Comparison ignores:
- Timing information (execution time, clock time)
- Step numbers
- CPU usage statistics
- Comments and metadata

## Test Requirements

Each test must:
1. Complete successfully (no errors)
2. Produce expected output files
3. Match gold files within specified tolerance (default: 1e-10)

## Tolerance Levels

- **1e-10** (default): Tight tolerance for most physics
- **1e-12**: Very strict, for analytical solutions
- **1e-8**: Relaxed, for cases with known minor differences

## Continuous Integration

For automated testing:
```bash
# Build
cd build_hh
cmake --build . --target fds_hh -j$(nproc)

# Run tests
cd ../test_cases
./run_tests.py

# Check exit code
if [ $? -eq 0 ]; then
    echo "All tests passed"
else
    echo "Some tests failed"
    exit 1
fi
```

## Troubleshooting

**Test fails with "FDS_HH not found":**
- Build fds_hh: `cd build_hh && cmake --build . --target fds_hh`

**Gold file missing:**
- Generate gold files: `./run_tests.py --generate-gold --test test_name`

**Comparison fails:**
- Check verbose output: `./run_tests.py --test test_name --verbose`
- Verify inputs match: `diff inputs/test.fds gold/test_name/test.fds`
- Increase tolerance if minor FP differences: `--tolerance 1e-8`

**Original FDS not available:**
- Use hedgehog for gold: `--use-hedgehog-gold`
- Or build original FDS from master branch

## File Cleanup

The test system automatically:
- Creates work directories in `run/`
- Cleans up intermediate files from `gold/` during generation
- Preserves only CSV output files in `gold/`

Manual cleanup:
```bash
# Remove all test run artifacts
rm -rf run/*

# Remove gold files (be careful!)
rm -rf gold/*
```

## Notes

- Gold files should be generated from the **master branch** using original FDS
- If comparing with hedgehog-generated gold, use `--use-hedgehog-gold` flag
- All tests use MPI with `--oversubscribe` to allow running on any system
- Test cases are designed to run quickly (< 5 seconds each)
