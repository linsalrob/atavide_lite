#!/usr/bin/env bash
# common.sh - Shared helper functions for atavide_lite pipeline
#
# Source this file in your scripts:
#   source "${ATAVIDE_ROOT}/lib/common.sh"
#   or
#   source "$(dirname "$0")/../lib/common.sh"

# Ensure this is sourced, not executed
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "Error: common.sh should be sourced, not executed directly" >&2
    exit 1
fi

# ============================================================================
# Logging and Error Handling
# ============================================================================

# Print error message to stderr and exit with non-zero status
# Usage: die "Error message"
die() {
    echo "ERROR: $*" >&2
    exit 1
}

# Log message with timestamp
# Usage: log "Processing sample X"
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*" >&2
}

# Log error message with timestamp (but don't exit)
# Usage: log_error "Warning: file not found"
log_error() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
}

# Log warning message with timestamp
# Usage: log_warn "This may cause issues"
log_warn() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] WARNING: $*" >&2
}

# ============================================================================
# Validation Functions
# ============================================================================

# Check if a command exists in PATH
# Usage: require_cmd fastp minimap2 samtools
require_cmd() {
    local missing=0
    for cmd in "$@"; do
        if ! command -v "$cmd" &> /dev/null; then
            log_error "Required command not found: $cmd"
            missing=1
        fi
    done
    if [[ $missing -eq 1 ]]; then
        die "Missing required commands. Please install them or load appropriate modules."
    fi
}

# Check if a file exists and is readable
# Usage: require_file "$input_file" "$reference_genome"
require_file() {
    local missing=0
    for file in "$@"; do
        if [[ ! -f "$file" ]]; then
            log_error "Required file not found: $file"
            missing=1
        elif [[ ! -r "$file" ]]; then
            log_error "Required file not readable: $file"
            missing=1
        fi
    done
    if [[ $missing -eq 1 ]]; then
        die "Missing or unreadable required files."
    fi
}

# Check if a directory exists, optionally create it
# Usage: require_dir "$output_dir" [create]
require_dir() {
    local dir="$1"
    local create="${2:-}"
    
    if [[ ! -d "$dir" ]]; then
        if [[ "$create" == "create" ]]; then
            log "Creating directory: $dir"
            mkdir -p "$dir" || die "Failed to create directory: $dir"
        else
            die "Required directory not found: $dir"
        fi
    fi
}

# Check if output file exists and has non-zero size
# Usage: check_nonempty "$output_file" || die "Output is empty"
check_nonempty() {
    local file="$1"
    if [[ ! -f "$file" ]]; then
        return 1
    fi
    if [[ ! -s "$file" ]]; then
        return 1
    fi
    return 0
}

# Check if a variable is set and non-empty
# Usage: require_var SAMPLE_NAME HOSTFILE
require_var() {
    local missing=0
    for var in "$@"; do
        if [[ -z "${!var:-}" ]]; then
            log_error "Required variable not set: $var"
            missing=1
        fi
    done
    if [[ $missing -eq 1 ]]; then
        die "Missing required environment variables."
    fi
}

# ============================================================================
# Configuration Loading
# ============================================================================

# Load configuration from paths.env file
# Usage: load_config [path_to_config]
load_config() {
    local config_file="${1:-config/paths.env}"
    
    # Try multiple possible locations
    local search_paths=(
        "$config_file"
        "$(pwd)/config/paths.env"
        "$(pwd)/../config/paths.env"
        "${ATAVIDE_ROOT:-}/config/paths.env"
    )
    
    local found=0
    for path in "${search_paths[@]}"; do
        if [[ -f "$path" ]]; then
            log "Loading configuration from: $path"
            # shellcheck disable=SC1090
            source "$path" || die "Failed to source config file: $path"
            found=1
            break
        fi
    done
    
    if [[ $found -eq 0 ]]; then
        log_error "Configuration file not found."
        log_error "Searched locations:"
        for path in "${search_paths[@]}"; do
            log_error "  - $path"
        done
        log_error ""
        log_error "Please create config/paths.env based on config/paths.env.example:"
        log_error "  cp config/paths.env.example config/paths.env"
        log_error "  # Edit config/paths.env with your paths"
        die "Configuration file not found"
    fi
    
    # Verify config was loaded
    if [[ "${ATAVIDE_CONFIG_LOADED:-0}" != "1" ]]; then
        log_warn "Config file loaded but ATAVIDE_CONFIG_LOADED not set"
    fi
}

# ============================================================================
# File Processing Helpers
# ============================================================================

# Count reads in a fastq.gz file
# Usage: read_count=$(count_fastq_reads "$file")
count_fastq_reads() {
    local file="$1"
    require_file "$file"
    
    if [[ "$file" == *.gz ]]; then
        zcat "$file" | echo $(($(wc -l)/4))
    else
        echo $(($(wc -l < "$file")/4))
    fi
}

# Count sequences in a fasta file (plain or gzipped)
# Usage: seq_count=$(count_fasta_sequences "$file")
count_fasta_sequences() {
    local file="$1"
    require_file "$file"
    
    if [[ "$file" == *.gz ]]; then
        zgrep -c '^>' "$file"
    else
        grep -c '^>' "$file"
    fi
}

# Get file size in human-readable format
# Usage: size=$(file_size "$file")
file_size() {
    local file="$1"
    require_file "$file"
    
    if command -v numfmt &> /dev/null; then
        numfmt --to=iec-i --suffix=B "$(stat -c%s "$file")"
    else
        ls -lh "$file" | awk '{print $5}'
    fi
}

# ============================================================================
# Slurm/PBS Helpers
# ============================================================================

# Get array task ID (works with both Slurm and PBS)
# Usage: task_id=$(get_array_task_id)
get_array_task_id() {
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        echo "${SLURM_ARRAY_TASK_ID}"
    elif [[ -n "${PBS_ARRAYID:-}" ]]; then
        echo "${PBS_ARRAYID}"
    else
        die "Not running as an array job (no SLURM_ARRAY_TASK_ID or PBS_ARRAYID)"
    fi
}

# Get job ID (works with both Slurm and PBS)
# Usage: job_id=$(get_job_id)
get_job_id() {
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
        echo "${SLURM_JOB_ID}"
    elif [[ -n "${PBS_JOBID:-}" ]]; then
        echo "${PBS_JOBID}"
    else
        echo "unknown"
    fi
}

# Detect scheduler type
# Usage: scheduler=$(detect_scheduler)
detect_scheduler() {
    if command -v sbatch &> /dev/null; then
        echo "slurm"
    elif command -v qsub &> /dev/null; then
        echo "pbs"
    else
        echo "unknown"
    fi
}

# ============================================================================
# System Detection
# ============================================================================

# Detect fast local storage (HPC system-specific)
# Usage: fast_storage=$(get_fast_storage)
get_fast_storage() {
    if [[ -n "${PBS_JOBFS:-}" && -d "${PBS_JOBFS}" ]]; then
        # NCI Gadi
        echo "${PBS_JOBFS}"
    elif [[ -n "${BGFS:-}" && -d "${BGFS}" ]]; then
        # Flinders Deepthought
        echo "${BGFS}"
    elif [[ -n "${SCRATCH_ROOT:-}" ]]; then
        # From config
        echo "${SCRATCH_ROOT}/tmp"
    elif [[ -d "/scratch" ]]; then
        # Pawsey or similar
        echo "/scratch/${PROJECT_NAME:-}/${USER}/tmp"
    else
        # Fallback
        echo "${TMPDIR:-/tmp}"
    fi
}

# ============================================================================
# Conda/Module Helpers
# ============================================================================

# Activate conda environment safely
# Usage: activate_conda "$ATAVIDE_CONDA"
activate_conda() {
    local env_name="$1"
    
    if [[ -z "$env_name" ]]; then
        die "Conda environment name required"
    fi
    
    # Initialize conda for bash
    if [[ -f "${CONDA_EXE%/bin/conda}/etc/profile.d/conda.sh" ]]; then
        # shellcheck disable=SC1091
        source "${CONDA_EXE%/bin/conda}/etc/profile.d/conda.sh"
    elif command -v conda &> /dev/null; then
        eval "$(conda shell.bash hook)"
    else
        die "Conda not found. Please install conda or set CONDA_EXE"
    fi
    
    log "Activating conda environment: $env_name"
    conda activate "$env_name" || die "Failed to activate conda environment: $env_name"
}

# ============================================================================
# Output Summary Helpers
# ============================================================================

# Print a summary box
# Usage: print_summary "Title" "Line 1" "Line 2" ...
print_summary() {
    local title="$1"
    shift
    
    echo "============================================================" >&2
    echo "  $title" >&2
    echo "============================================================" >&2
    for line in "$@"; do
        echo "  $line" >&2
    done
    echo "============================================================" >&2
}

# Print script completion message
# Usage: print_completion "$start_time"
print_completion() {
    local start_time="$1"
    local end_time
    end_time="$(date +%s)"
    local duration=$((end_time - start_time))
    
    log "Completed successfully in ${duration} seconds"
}

# ============================================================================
# Initialization
# ============================================================================

# Common initialization for all scripts
# Usage: init_script (call at start of your script)
init_script() {
    # Record start time
    export SCRIPT_START_TIME
    SCRIPT_START_TIME="$(date +%s)"
    
    log "=========================================="
    log "Script: ${0}"
    log "Host: $(hostname)"
    log "User: ${USER}"
    log "Date: $(date)"
    log "Job ID: $(get_job_id)"
    log "=========================================="
}

# Export functions for subshells if needed
export -f die log log_error log_warn
export -f require_cmd require_file require_dir check_nonempty require_var
export -f count_fastq_reads count_fasta_sequences file_size
export -f get_array_task_id get_job_id detect_scheduler get_fast_storage
export -f print_summary print_completion

# Mark that common.sh has been loaded
export ATAVIDE_COMMON_LOADED=1

log "Loaded common.sh helper functions"
