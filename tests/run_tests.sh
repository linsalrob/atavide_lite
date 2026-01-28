#!/usr/bin/env bash
# Main test runner for atavide_lite test suite
# 
# Usage: ./run_tests.sh [--help] [--verbose]
# 
# This script runs all test suites for the atavide_lite pipeline

set -euo pipefail

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Options
VERBOSE=0

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --help|-h)
            echo "Usage: $0 [--help] [--verbose]"
            echo ""
            echo "Run all test suites for atavide_lite pipeline"
            echo ""
            echo "Options:"
            echo "  --help, -h     Show this help message"
            echo "  --verbose, -v  Show verbose output"
            exit 0
            ;;
        --verbose|-v)
            VERBOSE=1
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Track overall results
TOTAL_SUITES=0
PASSED_SUITES=0
FAILED_SUITES=0

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}atavide_lite Test Suite Runner${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Function to run a test suite
run_test_suite() {
    local suite_name="$1"
    local suite_path="$2"
    
    TOTAL_SUITES=$((TOTAL_SUITES + 1))
    
    echo -e "${YELLOW}Running: ${suite_name}${NC}"
    echo "---"
    
    if [[ $VERBOSE -eq 1 ]]; then
        if "$suite_path"; then
            PASSED_SUITES=$((PASSED_SUITES + 1))
            echo ""
        else
            FAILED_SUITES=$((FAILED_SUITES + 1))
            echo ""
        fi
    else
        if "$suite_path" > /tmp/test_output.log 2>&1; then
            PASSED_SUITES=$((PASSED_SUITES + 1))
            echo -e "${GREEN}✓ ${suite_name} passed${NC}"
            echo ""
        else
            FAILED_SUITES=$((FAILED_SUITES + 1))
            echo -e "${RED}✗ ${suite_name} failed${NC}"
            echo "Output:"
            cat /tmp/test_output.log
            echo ""
        fi
    fi
}

# Run test suites
cd "${REPO_ROOT}"

# Test 1: lib/common.sh
if [[ -f "${SCRIPT_DIR}/test_common.sh" ]]; then
    chmod +x "${SCRIPT_DIR}/test_common.sh"
    run_test_suite "Bash Helper Functions (lib/common.sh)" "${SCRIPT_DIR}/test_common.sh"
fi

# Test 2: vamb_create_fasta_clusters.py
if [[ -f "${SCRIPT_DIR}/test_vamb_create_fasta_clusters.py" ]]; then
    chmod +x "${SCRIPT_DIR}/test_vamb_create_fasta_clusters.py"
    run_test_suite "VAMB Create FASTA Clusters" "${SCRIPT_DIR}/test_vamb_create_fasta_clusters.py"
fi

# Test 3: Configuration file validation
if [[ -f "${SCRIPT_DIR}/test_config_files.sh" ]]; then
    chmod +x "${SCRIPT_DIR}/test_config_files.sh"
    run_test_suite "Configuration Files" "${SCRIPT_DIR}/test_config_files.sh"
fi

# Print summary
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Test Suite Summary${NC}"
echo -e "${BLUE}========================================${NC}"
echo "Total test suites: ${TOTAL_SUITES}"
echo -e "${GREEN}Passed: ${PASSED_SUITES}${NC}"

if [[ $FAILED_SUITES -gt 0 ]]; then
    echo -e "${RED}Failed: ${FAILED_SUITES}${NC}"
    echo ""
    echo -e "${RED}Some test suites failed. Please review the output above.${NC}"
    exit 1
else
    echo ""
    echo -e "${GREEN}All test suites passed!${NC}"
    exit 0
fi
