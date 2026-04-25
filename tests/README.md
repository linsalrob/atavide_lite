# Test Suite for atavide_lite

This directory contains the test suite for the atavide_lite pipeline improvements.

## Overview

The test suite validates:
1. **lib/common.sh** - Bash helper functions
2. **bin/vamb_create_fasta_clusters.py** - VAMB binning script
3. **Configuration files** - paths.env.example and samples.tsv.example

## Quick Start

Run all tests:
```bash
cd tests
./run_tests.sh
```

Run with verbose output:
```bash
./run_tests.sh --verbose
```

Run individual test suites:
```bash
./test_common.sh                      # Test Bash helpers
./test_vamb_create_fasta_clusters.py  # Test VAMB script
./test_config_files.sh                # Test config files
```

## Test Suites

### 1. test_common.sh

Tests the Bash helper functions in `lib/common.sh`:

**Logging Functions**
- `log()` - Timestamped logging
- `log_error()` - Error logging
- `log_warn()` - Warning logging
- `die()` - Fatal error with exit

**Validation Functions**
- `require_cmd()` - Check command exists
- `require_file()` - Check file exists and readable
- `require_dir()` - Check/create directory
- `check_nonempty()` - Verify non-empty file
- `require_var()` - Check environment variable set

**File Processing Functions**
- `count_fastq_reads()` - Count FASTQ reads
- `count_fasta_sequences()` - Count FASTA sequences
- `file_size()` - Get human-readable size

**HPC Helper Functions**
- `get_array_task_id()` - Get Slurm/PBS array ID
- `get_job_id()` - Get job ID
- `detect_scheduler()` - Detect Slurm/PBS/unknown
- `get_fast_storage()` - Detect fast local storage

**Configuration Functions**
- `load_config()` - Load config/paths.env

**Output Functions**
- `print_summary()` - Formatted summary box
- `print_completion()` - Completion with duration
- `init_script()` - Initialize with logging

**Coverage**: ~30 test cases

### 2. test_vamb_create_fasta_clusters.py

Tests the VAMB cluster FASTA creation script:

**Basic Functionality**
- Help message display
- Argument validation
- Basic cluster FASTA creation
- Compressed output (gzip)

**Advanced Features**
- Minimum size filtering (`-m` parameter)
- Verbose output (`-v` parameter)
- Multi-cluster handling

**Error Handling**
- Missing input files
- Malformed cluster files
- Invalid arguments

**Coverage**: ~8 test cases

### 3. test_config_files.sh

Tests the configuration file templates:

**paths.env.example**
- File exists
- Contains required variables
- Valid Bash syntax
- Sourceable without errors

**samples.tsv.example**
- File exists
- Has proper header columns
- Contains paired-end examples
- Contains single-end examples
- Uses tab delimiters

**Coverage**: ~15 test cases

## Test Output

Each test suite produces colored output:
- ✓ **Green**: Test passed
- ✗ **Red**: Test failed
- ⚠ **Yellow**: Warning

Example output:
```
Testing: Logging Functions
---
✓ log() produces output with message
✓ log() includes timestamp
✓ log_error() produces error output
✓ log_warn() produces warning output

========================================
Test Results
========================================
Total tests run: 30
Passed: 30
All tests passed!
```

## Requirements

### Bash Tests
- Bash 4.0+
- Standard Unix utilities (grep, awk, etc.)

### Python Tests
- Python 3.9+
- No external dependencies (uses stdlib only)
- atavide_lib module (for full integration)

## Adding New Tests

### Adding a Bash Test

1. Create test function:
```bash
test_my_function() {
    # Test logic
    some_command
}
assert_success "Test description" test_my_function
```

2. Use assert helpers:
- `assert_success "desc" command` - Command should succeed
- `assert_failure "desc" command` - Command should fail
- `assert_equals "desc" expected actual` - Values should match
- `assert_contains "desc" substring text` - Text contains substring

### Adding a Python Test

1. Add method to TestVambCreateFastaClusters class:
```python
def test_my_feature(self):
    """Test description"""
    tmpdir = self.setup()
    try:
        # Test logic
        result = self.run_script(...)
        assert_equals("desc", expected, actual)
    finally:
        self.teardown()
```

2. Add to `run_all_tests()`:
```python
self.test_my_feature()
```

### Adding a New Test Suite

1. Create test script: `test_newfeature.sh` or `test_newfeature.py`
2. Add to `run_tests.sh`:
```bash
if [[ -f "${SCRIPT_DIR}/test_newfeature.sh" ]]; then
    chmod +x "${SCRIPT_DIR}/test_newfeature.sh"
    run_test_suite "New Feature" "${SCRIPT_DIR}/test_newfeature.sh"
fi
```

## Continuous Integration

These tests can be integrated into GitHub Actions:

```yaml
name: Test Suite

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Run tests
        run: |
          cd tests
          ./run_tests.sh
```

## Test Coverage Goals

- **lib/common.sh**: 100% function coverage ✅
- **vamb_create_fasta_clusters.py**: Core functionality ✅
- **Configuration files**: Syntax and structure ✅
- **Integration tests**: Example scripts (future)

## Known Limitations

1. **No VAMB module**: Tests don't require VAMB installed (by design)
2. **No HPC environment**: Some HPC-specific features tested with mocks
3. **No real data**: Tests use synthetic small datasets
4. **No integration tests**: Full pipeline integration testing is manual

## Troubleshooting

**Tests fail with "command not found"**
- Ensure you're running from the tests/ directory
- Check that scripts have execute permissions: `chmod +x test_*.sh`

**Python tests fail with import errors**
- Ensure atavide_lib is in the repository root
- Check Python version: `python3 --version` (needs 3.9+)

**Bash tests fail on macOS**
- Some commands may differ (e.g., `stat` syntax)
- Tests are primarily designed for Linux

## Future Enhancements

Potential additions to the test suite:
- Integration tests for full pipeline runs
- Performance/benchmark tests
- Database compatibility tests
- Cluster script validation (Slurm/PBS syntax)
- Documentation link checking
- Shellcheck integration for all scripts

## Contributing

When adding new features to atavide_lite:
1. Add corresponding tests to the test suite
2. Run tests locally before committing
3. Ensure all tests pass
4. Update this README if adding new test suites

## License

Same as the main atavide_lite project (MIT).
