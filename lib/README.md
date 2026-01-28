# lib/ - Shared Helper Functions

This directory contains shared Bash functions for use across atavide_lite scripts.

## Files

### common.sh

Shared helper functions for pipeline scripts. See the full file for detailed documentation of each function.

**Usage in scripts:**
```bash
#!/usr/bin/env bash
set -euo pipefail

# Source common.sh
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../lib/common.sh" || {
    echo "ERROR: Could not source lib/common.sh" >&2
    exit 1
}

# Initialize logging
init_script

# Use helper functions
require_cmd fastp minimap2 samtools
require_file "input.fastq.gz"
require_dir "output" create

log "Processing started"
# ... your code here ...
print_completion "${SCRIPT_START_TIME}"
```

## Available Functions

### Logging and Error Handling
- `die "message"` - Print error and exit
- `log "message"` - Log with timestamp
- `log_error "message"` - Log error (doesn't exit)
- `log_warn "message"` - Log warning

### Validation
- `require_cmd cmd1 cmd2 ...` - Check commands exist in PATH
- `require_file file1 file2 ...` - Check files exist and are readable
- `require_dir dir [create]` - Check directory exists, optionally create
- `require_var VAR1 VAR2 ...` - Check variables are set
- `check_nonempty file` - Check file exists and has size > 0

### Configuration
- `load_config [path]` - Load config/paths.env file

### File Processing
- `count_fastq_reads file.fastq.gz` - Count reads in fastq file
- `count_fasta_sequences file.fasta.gz` - Count sequences in fasta
- `file_size file` - Get human-readable file size

### HPC Helpers
- `get_array_task_id` - Get Slurm/PBS array task ID
- `get_job_id` - Get Slurm/PBS job ID
- `detect_scheduler` - Detect slurm/pbs/unknown
- `get_fast_storage` - Detect fast local storage path

### Conda
- `activate_conda env_name` - Safely activate conda environment

### Output
- `print_summary "Title" "Line1" "Line2" ...` - Print formatted summary box
- `print_completion start_time` - Print completion message with duration
- `init_script` - Initialize script with logging (sets SCRIPT_START_TIME)

## Example: Enhanced Script

See `pawsey_shortread/fastp_enhanced.slurm` for a complete example of using lib/common.sh in a cluster script.

Key improvements when using lib/common.sh:
1. **Better error handling** - Fail fast with clear messages
2. **Validation** - Check prerequisites before running
3. **Logging** - Timestamped messages for debugging
4. **Portability** - Works across different HPC systems
5. **Consistency** - Same patterns across all scripts

## Migration Guide

To update an existing script to use lib/common.sh:

1. **Add source statement** at top (after shebang and SBATCH headers):
   ```bash
   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
   source "${SCRIPT_DIR}/../lib/common.sh"
   ```

2. **Replace manual checks** with helper functions:
   ```bash
   # Before:
   if [[ ! -f "$file" ]]; then
       echo "Error: file not found" >&2
       exit 1
   fi
   
   # After:
   require_file "$file"
   ```

3. **Add initialization** at script start:
   ```bash
   init_script
   ```

4. **Replace echo with log**:
   ```bash
   # Before:
   echo "Processing sample X" >&2
   
   # After:
   log "Processing sample X"
   ```

5. **Add completion message** at end:
   ```bash
   print_completion "${SCRIPT_START_TIME}"
   ```

## Design Principles

- **POSIX-compatible where possible** - Works on any Bash 4.0+
- **No external dependencies** - Only uses standard Unix tools
- **Fail-fast** - Validate early, fail with clear messages
- **Logging** - Always log to stderr, leave stdout for data
- **Portable** - Detects HPC system differences automatically

## Testing

Test that common.sh loads correctly:
```bash
source lib/common.sh
log "Test message"
require_cmd bash
echo "Success!"
```

## Contributing

When adding new functions:
1. Document with comments above the function
2. Include usage example in comment
3. Use consistent naming (verb_noun pattern)
4. Export functions if needed for subshells
5. Update this README
