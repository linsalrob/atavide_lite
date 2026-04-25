#!/usr/bin/env bash
# Test suite for configuration files
# 
# Usage: ./test_config_files.sh
# 
# This script validates the example configuration files

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

# Test counters
TESTS_RUN=0
TESTS_PASSED=0
TESTS_FAILED=0

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Helper functions
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

assert_file_exists() {
    TESTS_RUN=$((TESTS_RUN + 1))
    local test_name="$1"
    local filepath="$2"
    
    if [[ -f "$filepath" ]]; then
        echo -e "${GREEN}✓${NC} $test_name"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        return 0
    else
        echo -e "${RED}✗${NC} $test_name (file not found: $filepath)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
}

assert_file_contains() {
    TESTS_RUN=$((TESTS_RUN + 1))
    local test_name="$1"
    local filepath="$2"
    local pattern="$3"
    
    if grep -q "$pattern" "$filepath" 2>/dev/null; then
        echo -e "${GREEN}✓${NC} $test_name"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        return 0
    else
        echo -e "${RED}✗${NC} $test_name"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
}

# ============================================================================
# Test: paths.env.example
# ============================================================================
echo "Testing: config/paths.env.example"
echo "---"

PATHS_ENV="${REPO_ROOT}/config/paths.env.example"

assert_file_exists "paths.env.example exists" "$PATHS_ENV"

# Check for required variables
assert_file_contains "Contains SCRATCH_ROOT" "$PATHS_ENV" "SCRATCH_ROOT"
assert_file_contains "Contains TMPDIR" "$PATHS_ENV" "TMPDIR"
assert_file_contains "Contains FASTQ_DIR" "$PATHS_ENV" "FASTQ_DIR"
assert_file_contains "Contains HOSTFILE" "$PATHS_ENV" "HOSTFILE"
assert_file_contains "Contains MMSEQS_DB_DIR" "$PATHS_ENV" "MMSEQS_DB_DIR"
assert_file_contains "Contains THREADS" "$PATHS_ENV" "THREADS"
assert_file_contains "Contains ATAVIDE_CONDA" "$PATHS_ENV" "ATAVIDE_CONDA"
assert_file_contains "Contains ATAVIDE_CONFIG_LOADED marker" "$PATHS_ENV" "ATAVIDE_CONFIG_LOADED"

# Check that it's sourceable (syntax check)
test_paths_env_sourceable() {
    bash -n "$PATHS_ENV"
}
assert_success "paths.env.example has valid bash syntax" test_paths_env_sourceable

# ============================================================================
# Test: samples.tsv.example
# ============================================================================
echo ""
echo "Testing: config/samples.tsv.example"
echo "---"

SAMPLES_TSV="${REPO_ROOT}/config/samples.tsv.example"

assert_file_exists "samples.tsv.example exists" "$SAMPLES_TSV"

# Check for header
assert_file_contains "Contains header with sample_id" "$SAMPLES_TSV" "sample_id"
assert_file_contains "Contains header with r1" "$SAMPLES_TSV" "r1"
assert_file_contains "Contains header with r2" "$SAMPLES_TSV" "r2"

# Check for example rows
assert_file_contains "Contains paired-end example" "$SAMPLES_TSV" "sample001"
assert_file_contains "Contains single-end example" "$SAMPLES_TSV" "sample_nanopore"

# Check it's tab-delimited
test_tsv_format() {
    grep -q $'\t' "$SAMPLES_TSV"
}
assert_success "samples.tsv.example uses tab delimiters" test_tsv_format

# ============================================================================
# Test: Config directory structure
# ============================================================================
echo ""
echo "Testing: Configuration directory structure"
echo "---"

# Check config directory exists
test_config_dir_exists() {
    [[ -d "${REPO_ROOT}/config" ]]
}
assert_success "config directory exists" test_config_dir_exists

# Check that README or guidance exists
test_config_has_docs() {
    [[ -f "${REPO_ROOT}/config/paths.env.example" ]] && \
    [[ -f "${REPO_ROOT}/config/samples.tsv.example" ]]
}
assert_success "Config directory has example files" test_config_has_docs

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
