#!/usr/bin/env bash
# Test suite for lib/common.sh
# 
# Usage: ./test_common.sh
# 
# This script tests the helper functions in lib/common.sh

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test counters
TESTS_RUN=0
TESTS_PASSED=0
TESTS_FAILED=0

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Create temporary test directory
TEST_TMPDIR=$(mktemp -d -t atavide_test.XXXXXX)
trap 'rm -rf "${TEST_TMPDIR}"' EXIT

# Source common.sh in a controlled way
export ATAVIDE_ROOT="${REPO_ROOT}"
cd "${TEST_TMPDIR}"

# Redirect sourcing errors to not pollute test output
if ! source "${REPO_ROOT}/lib/common.sh" 2>/dev/null; then
    echo -e "${RED}✗ FATAL: Could not source lib/common.sh${NC}"
    exit 1
fi

# Helper functions for testing
assert_success() {
    TESTS_RUN=$((TESTS_RUN + 1))
    local test_name="$1"
    shift
    
    if "$@" &>/dev/null; then
        echo -e "${GREEN}✓${NC} $test_name"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        return 0
    else
        echo -e "${RED}✗${NC} $test_name"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
}

assert_failure() {
    TESTS_RUN=$((TESTS_RUN + 1))
    local test_name="$1"
    shift
    
    if ! "$@" &>/dev/null; then
        echo -e "${GREEN}✓${NC} $test_name"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        return 0
    else
        echo -e "${RED}✗${NC} $test_name"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
}

assert_equals() {
    TESTS_RUN=$((TESTS_RUN + 1))
    local test_name="$1"
    local expected="$2"
    local actual="$3"
    
    if [[ "$expected" == "$actual" ]]; then
        echo -e "${GREEN}✓${NC} $test_name"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        return 0
    else
        echo -e "${RED}✗${NC} $test_name (expected: '$expected', got: '$actual')"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
}

assert_contains() {
    TESTS_RUN=$((TESTS_RUN + 1))
    local test_name="$1"
    local substring="$2"
    local text="$3"
    
    if [[ "$text" == *"$substring"* ]]; then
        echo -e "${GREEN}✓${NC} $test_name"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        return 0
    else
        echo -e "${RED}✗${NC} $test_name (substring '$substring' not found)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
}

# ============================================================================
# Test: Logging Functions
# ============================================================================
echo "Testing: Logging Functions"
echo "---"

test_log_output() {
    local output
    output=$(log "test message" 2>&1)
    [[ "$output" == *"test message"* ]]
}
assert_success "log() produces output with message" test_log_output

test_log_timestamp() {
    local output
    output=$(log "test" 2>&1)
    # Check for timestamp pattern [YYYY-MM-DD HH:MM:SS]
    [[ "$output" =~ \[[0-9]{4}-[0-9]{2}-[0-9]{2}\ [0-9]{2}:[0-9]{2}:[0-9]{2}\] ]]
}
assert_success "log() includes timestamp" test_log_timestamp

test_log_error_output() {
    local output
    output=$(log_error "error message" 2>&1)
    [[ "$output" == *"ERROR: error message"* ]]
}
assert_success "log_error() produces error output" test_log_error_output

test_log_warn_output() {
    local output
    output=$(log_warn "warning message" 2>&1)
    [[ "$output" == *"WARNING: warning message"* ]]
}
assert_success "log_warn() produces warning output" test_log_warn_output

# ============================================================================
# Test: Validation Functions
# ============================================================================
echo ""
echo "Testing: Validation Functions"
echo "---"

test_require_cmd_success() {
    require_cmd bash
}
assert_success "require_cmd() succeeds for existing command" test_require_cmd_success

test_require_cmd_failure() {
    (require_cmd nonexistent_command_xyz_123)
}
assert_failure "require_cmd() fails for missing command" test_require_cmd_failure

test_require_file_success() {
    echo "test" > "${TEST_TMPDIR}/testfile.txt"
    require_file "${TEST_TMPDIR}/testfile.txt"
}
assert_success "require_file() succeeds for existing file" test_require_file_success

test_require_file_failure() {
    (require_file "${TEST_TMPDIR}/nonexistent.txt")
}
assert_failure "require_file() fails for missing file" test_require_file_failure

test_require_dir_existing() {
    mkdir -p "${TEST_TMPDIR}/testdir"
    require_dir "${TEST_TMPDIR}/testdir"
}
assert_success "require_dir() succeeds for existing directory" test_require_dir_existing

test_require_dir_create() {
    require_dir "${TEST_TMPDIR}/newdir" create
    [[ -d "${TEST_TMPDIR}/newdir" ]]
}
assert_success "require_dir() creates directory with 'create' option" test_require_dir_create

test_check_nonempty_success() {
    echo "content" > "${TEST_TMPDIR}/nonempty.txt"
    check_nonempty "${TEST_TMPDIR}/nonempty.txt"
}
assert_success "check_nonempty() succeeds for non-empty file" test_check_nonempty_success

test_check_nonempty_empty() {
    touch "${TEST_TMPDIR}/empty.txt"
    check_nonempty "${TEST_TMPDIR}/empty.txt"
}
assert_failure "check_nonempty() fails for empty file" test_check_nonempty_empty

test_require_var_success() {
    export TEST_VAR="value"
    require_var TEST_VAR
}
assert_success "require_var() succeeds for set variable" test_require_var_success

test_require_var_failure() {
    unset TEST_VAR_UNSET 2>/dev/null || true
    (require_var TEST_VAR_UNSET)
}
assert_failure "require_var() fails for unset variable" test_require_var_failure

# ============================================================================
# Test: File Processing Functions
# ============================================================================
echo ""
echo "Testing: File Processing Functions"
echo "---"

test_count_fastq_reads() {
    # Create a test fastq file (4 lines = 1 read)
    cat > "${TEST_TMPDIR}/test.fastq" <<EOF
@read1
ACGT
+
IIII
@read2
TGCA
+
IIII
EOF
    local count
    count=$(count_fastq_reads "${TEST_TMPDIR}/test.fastq")
    [[ "$count" == "2" ]]
}
assert_success "count_fastq_reads() counts reads correctly" test_count_fastq_reads

test_count_fasta_sequences() {
    cat > "${TEST_TMPDIR}/test.fasta" <<EOF
>seq1
ACGT
>seq2
TGCA
>seq3
AAAA
EOF
    local count
    count=$(count_fasta_sequences "${TEST_TMPDIR}/test.fasta")
    [[ "$count" == "3" ]]
}
assert_success "count_fasta_sequences() counts sequences correctly" test_count_fasta_sequences

test_file_size() {
    echo "test content" > "${TEST_TMPDIR}/sizefile.txt"
    local size
    size=$(file_size "${TEST_TMPDIR}/sizefile.txt")
    [[ -n "$size" ]]  # Just check it returns something
}
assert_success "file_size() returns file size" test_file_size

# ============================================================================
# Test: HPC Helper Functions
# ============================================================================
echo ""
echo "Testing: HPC Helper Functions"
echo "---"

test_detect_scheduler() {
    local scheduler
    scheduler=$(detect_scheduler)
    [[ -n "$scheduler" ]]  # Returns something (slurm, pbs, or unknown)
}
assert_success "detect_scheduler() returns a value" test_detect_scheduler

test_get_fast_storage() {
    local storage
    storage=$(get_fast_storage)
    [[ -n "$storage" ]]  # Returns a path
}
assert_success "get_fast_storage() returns a path" test_get_fast_storage

# Test get_job_id (should return "unknown" when not in a job)
test_get_job_id() {
    local job_id
    job_id=$(get_job_id)
    [[ "$job_id" == "unknown" ]]
}
assert_success "get_job_id() returns 'unknown' outside of job" test_get_job_id

# ============================================================================
# Test: Output Functions
# ============================================================================
echo ""
echo "Testing: Output Functions"
echo "---"

test_print_summary() {
    local output
    output=$(print_summary "Test Title" "Line 1" "Line 2" 2>&1)
    [[ "$output" == *"Test Title"* ]] && [[ "$output" == *"Line 1"* ]]
}
assert_success "print_summary() prints title and lines" test_print_summary

test_init_script() {
    unset SCRIPT_START_TIME
    init_script >/dev/null 2>&1
    [[ -n "${SCRIPT_START_TIME:-}" ]]
}
assert_success "init_script() sets SCRIPT_START_TIME" test_init_script

test_print_completion() {
    SCRIPT_START_TIME=$(date +%s)
    sleep 1
    local output
    output=$(print_completion "${SCRIPT_START_TIME}" 2>&1)
    [[ "$output" == *"Completed successfully"* ]]
}
assert_success "print_completion() prints completion message" test_print_completion

# ============================================================================
# Test: Config Loading
# ============================================================================
echo ""
echo "Testing: Configuration Loading"
echo "---"

test_load_config_missing() {
    cd "${TEST_TMPDIR}"
    (load_config)
}
assert_failure "load_config() fails when config not found" test_load_config_missing

test_load_config_success() {
    mkdir -p "${TEST_TMPDIR}/config"
    cat > "${TEST_TMPDIR}/config/paths.env" <<EOF
export TEST_CONFIG_VAR="loaded"
export ATAVIDE_CONFIG_LOADED=1
EOF
    cd "${TEST_TMPDIR}"
    unset TEST_CONFIG_VAR 2>/dev/null || true
    load_config >/dev/null 2>&1
    [[ "${TEST_CONFIG_VAR:-}" == "loaded" ]]
}
assert_success "load_config() loads config file" test_load_config_success

# ============================================================================
# Summary
# ============================================================================
echo ""
echo "========================================"
echo "Test Results"
echo "========================================"
echo "Total tests run: $TESTS_RUN"
echo -e "${GREEN}Passed: $TESTS_PASSED${NC}"
if [[ $TESTS_FAILED -gt 0 ]]; then
    echo -e "${RED}Failed: $TESTS_FAILED${NC}"
    echo ""
    echo "Some tests failed. Please review the output above."
    exit 1
else
    echo -e "${GREEN}All tests passed!${NC}"
    exit 0
fi
