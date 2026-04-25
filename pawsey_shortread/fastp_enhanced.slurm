#!/usr/bin/env bash
#SBATCH --job-name=fastp
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH -o slurm_output/fastp_slurm/fastp-%A_%a.out
#SBATCH -e slurm_output/fastp_slurm/fastp-%A_%a.err

# fastp QC/Adapter Trimming Script
# Enhanced with lib/common.sh helper functions
# 
# This is an EXAMPLE of how to use lib/common.sh in cluster scripts.
# The original fastp.slurm is still available and functional.

set -euo pipefail

# Source common helper functions
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "${SCRIPT_DIR}/../lib/common.sh" || {
    echo "ERROR: Could not source lib/common.sh" >&2
    exit 1
}

# Initialize script with logging
init_script

# Check for required environment variables
if [[ ! -n "${ATAVIDE_CONDA:-}" ]]; then
    die "ATAVIDE_CONDA environment variable not set. Please define the location of your ATAVIDE_LITE conda installation"
fi

# Activate conda environment
log "Activating conda environment: ${ATAVIDE_CONDA}"
eval "$(conda shell.bash hook)"
conda activate "${ATAVIDE_CONDA}" || die "Failed to activate conda environment"

# Validate required files
require_file "R1_reads.txt" || {
    log_error "R1_reads.txt not found"
    log_error "Create it with: find fastq -name '*R1*' -printf \"%f\n\" > R1_reads.txt"
    exit 2
}

require_file "DEFINITIONS.sh" || {
    log_error "DEFINITIONS.sh not found"
    log_error "Please create DEFINITIONS.sh with SOURCE, FILEEND, HOSTREMOVED variables"
    exit 2
}

# Source definitions
log "Loading DEFINITIONS.sh"
# shellcheck disable=SC1091
source DEFINITIONS.sh

# Validate required variables from DEFINITIONS
require_var SOURCE FILEEND

# Get array task ID
TASK_ID=$(get_array_task_id)
log "Processing array task: ${TASK_ID}"

# Get R1 and R2 file names
R1=$(head -n "${TASK_ID}" R1_reads.txt | tail -n 1)
R2="${R1/_R1/_R2}"

log "Processing sample: ${R1} and ${R2}"

# Validate input files exist
require_file "${SOURCE}/${R1}" "${SOURCE}/${R2}"

# Create output directories
require_dir "fastq_fastp" create
require_dir "fastp_output" create
require_dir "slurm_output/fastp_slurm" create

# Set adapter file path
ADAPTER_FILE="${HOME}/atavide_lite/adapters/IlluminaAdapters.fa"
if [[ ! -f "${ADAPTER_FILE}" ]]; then
    log_warn "Adapter file not found at ${ADAPTER_FILE}, using default adapters"
    ADAPTER_FILE=""
fi

# Run fastp
log "Running fastp for ${R1}"
FASTP_CMD="fastp -n 1 -l 100 -i ${SOURCE}/${R1} -I ${SOURCE}/${R2} -o fastq_fastp/${R1} -O fastq_fastp/${R2} -j fastp_output/${R1}.json -h fastp_output/${R1}.html --thread 16"
if [[ -n "${ADAPTER_FILE}" ]]; then
    FASTP_CMD="${FASTP_CMD} --adapter_fasta ${ADAPTER_FILE}"
fi

log "Command: ${FASTP_CMD}"
${FASTP_CMD} || die "fastp failed"

# Validate outputs
log "Validating outputs"
check_nonempty "fastq_fastp/${R1}" || die "Output R1 file is empty or missing"
check_nonempty "fastq_fastp/${R2}" || die "Output R2 file is empty or missing"
check_nonempty "fastp_output/${R1}.json" || log_warn "JSON report is empty or missing"

# Print summary
INPUT_R1_SIZE=$(file_size "${SOURCE}/${R1}")
OUTPUT_R1_SIZE=$(file_size "fastq_fastp/${R1}")

print_summary "fastp Completed Successfully" \
    "Sample: ${R1}" \
    "Input size: ${INPUT_R1_SIZE}" \
    "Output size: ${OUTPUT_R1_SIZE}" \
    "Report: fastp_output/${R1}.html"

print_completion "${SCRIPT_START_TIME}"
