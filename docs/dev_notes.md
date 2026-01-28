# Development Notes

This document provides guidance for developers working on the atavide_lite pipeline.

---

## Code Quality Tools

### ShellCheck

**ShellCheck** is a static analysis tool for shell scripts that catches common bugs and style issues.

#### Installation
```bash
# Ubuntu/Debian
apt-get install shellcheck

# macOS
brew install shellcheck

# Conda
conda install -c conda-forge shellcheck
```

#### Usage
```bash
# Check a single script
shellcheck path/to/script.sh

# Check all scripts in a directory
find pawsey_shortread -name "*.sh" -o -name "*.slurm" | xargs shellcheck

# Exclude specific warnings (if needed)
shellcheck -e SC1090 script.sh  # Exclude "Can't follow non-constant source"
```

#### Common Issues to Watch For

1. **SC2086**: Unquoted variables
   ```bash
   # Bad
   echo $VAR
   
   # Good
   echo "${VAR}"
   ```

2. **SC2154**: Variable referenced but not assigned
   ```bash
   # Check that variables are defined or sourced
   if [[ ! -n "$ATAVIDE_CONDA" ]]; then
       echo "Error: ATAVIDE_CONDA not set" >&2
       exit 1
   fi
   ```

3. **SC2046**: Quote to prevent word splitting
   ```bash
   # Bad
   for file in $(find . -name "*.fastq"); do
   
   # Good
   find . -name "*.fastq" -print0 | while IFS= read -r -d '' file; do
   ```

4. **SC2181**: Check exit code directly
   ```bash
   # Bad
   command
   if [ $? -eq 0 ]; then
   
   # Good
   if command; then
   ```

#### Integration with Development Workflow

**Pre-commit check** (optional):
```bash
#!/bin/bash
# Save as .git/hooks/pre-commit

echo "Running shellcheck on modified scripts..."
git diff --cached --name-only --diff-filter=ACM | \
    grep -E '\.(sh|bash|slurm)$' | \
    xargs shellcheck || exit 1
```

---

## Shell Scripting Standards

### Script Header
```bash
#!/usr/bin/env bash
# Description: What this script does
# Usage: script.sh [options]
#
# Author: Name
# Last modified: YYYY-MM-DD

set -euo pipefail
```

### Error Handling
```bash
# Define at top of script
function die() {
    echo "ERROR: $*" >&2
    exit 1
}

# Use for fatal errors
[[ -f "$INPUT_FILE" ]] || die "Input file not found: $INPUT_FILE"

# Trap for debugging
trap 'echo "Error on line $LINENO: $BASH_COMMAND"' ERR
```

### Variable Naming
```bash
# Use UPPERCASE for exported/environment variables
export ATAVIDE_CONDA=/path/to/env

# Use lowercase for local variables
input_file="data.txt"

# Use descriptive names
sample_count=10  # Good
x=10             # Bad
```

### Quoting
```bash
# Always quote variables
echo "${var}"

# Quote command substitutions
files="$(find . -name '*.txt')"

# Use arrays for multiple values
files=( file1.txt file2.txt file3.txt )
for file in "${files[@]}"; do
    echo "$file"
done
```

### Conditionals
```bash
# Use [[ ]] for tests (bash-specific, more features)
if [[ -f "$file" && -r "$file" ]]; then
    # Good
fi

# Not [ ] (POSIX, fewer features)

# Prefer explicit comparison
if [[ "$var" == "value" ]]; then  # Good
if [ "$var" = "value" ]; then     # POSIX compatible
```

### Loops and File Processing
```bash
# For files with spaces/special characters
find . -name "*.fastq" -print0 | while IFS= read -r -d '' file; do
    process "$file"
done

# Simple iteration
for file in *.txt; do
    [[ -e "$file" ]] || continue  # Handle no matches
    echo "$file"
done
```

---

## Python Code Standards

### Style Guide
- Follow **PEP 8** (enforced by tools below)
- Maximum line length: 100 characters (compromise for readability)
- Use type hints for function signatures

### Tools

#### Black (code formatter)
```bash
# Install
pip install black

# Format a file
black script.py

# Format all Python files
find bin -name "*.py" | xargs black

# Check without modifying
black --check script.py
```

#### Flake8 (linter)
```bash
# Install
pip install flake8

# Run on a file
flake8 script.py

# Configure in setup.cfg or .flake8
[flake8]
max-line-length = 100
exclude = .git,__pycache__,build,dist
ignore = E203,W503
```

#### mypy (type checker)
```bash
# Install
pip install mypy

# Check types
mypy script.py

# Relax for gradual adoption
mypy --ignore-missing-imports script.py
```

### Example Python Script Template
```python
#!/usr/bin/env python3
"""
Brief description of what this script does.

Longer description if needed.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional


def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '-i', '--input',
        type=Path,
        required=True,
        help='Input file path'
    )
    parser.add_argument(
        '-o', '--output',
        type=Path,
        required=True,
        help='Output file path'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    # Process
    try:
        process_file(args.input, args.output, verbose=args.verbose)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


def process_file(input_path: Path, output_path: Path, verbose: bool = False) -> None:
    """Process the input file and write to output."""
    if verbose:
        print(f"Processing {input_path}...", file=sys.stderr)
    
    # Implementation here
    pass


if __name__ == '__main__':
    main()
```

---

## Testing

### Manual Testing Checklist

When modifying scripts:

1. **Syntax check**:
   ```bash
   bash -n script.sh  # Check syntax without running
   ```

2. **Dry run** (if supported):
   ```bash
   ./script.sh --dry-run
   ```

3. **Test with minimal data**:
   - Use small test files (1000 reads)
   - Verify outputs match expected format
   - Check error handling

4. **Test on target system**:
   - Run on actual HPC system
   - Verify module/conda environment
   - Check resource usage

### Python Unit Tests (future)

For critical helper scripts, consider adding tests:

```python
# test_vamb_helpers.py
import pytest
from vamb_create_fasta_clusters import parse_clusters


def test_parse_clusters():
    """Test cluster parsing."""
    # Implementation
    pass
```

---

## Git Workflow

### Commit Messages
```
<type>: <short summary>

<longer description if needed>

- Bullet points for details
- What changed and why
```

**Types**: fix, feat, docs, style, refactor, test, chore

**Examples**:
```
fix: Handle VAMB v4 API changes in vamb_create_fasta.py

docs: Add directory contract documentation

feat: Add lib/common.sh with shared helper functions
```

### Branching
- **main**: Stable, tested code
- **feature branches**: `feature/description` or `fix/bug-description`
- **Test before merging**: Run on small dataset

---

## Common Patterns to Avoid

### 1. Unquoted Variables
```bash
# BAD
file=$1
cat $file

# GOOD
file="$1"
cat "$file"
```

### 2. Using `ls` for Iteration
```bash
# BAD
for file in $(ls *.txt); do

# GOOD
for file in *.txt; do
    [[ -e "$file" ]] || continue
```

### 3. Not Checking Exit Codes
```bash
# BAD
command_that_might_fail
next_command

# GOOD
if ! command_that_might_fail; then
    echo "Command failed" >&2
    exit 1
fi
next_command
```

### 4. Hardcoded Paths
```bash
# BAD
DB=/scratch/pawsey0001/user/Databases/UniRef50

# GOOD
DB="${MMSEQS_DB_DIR}/${MMSEQS_DB_NAME}"
```

### 5. Missing Input Validation
```bash
# BAD
input_file=$1
cat "$input_file"

# GOOD
input_file="${1:?Error: Input file required}"
if [[ ! -f "$input_file" ]]; then
    echo "Error: File not found: $input_file" >&2
    exit 1
fi
cat "$input_file"
```

---

## Documentation Standards

### Code Comments
```bash
# Brief explanation of what the next block does
# Not obvious "what" but explaining "why"

# Bad comment:
# Loop through files
for file in *.txt; do

# Good comment:
# Process each sample separately to allow resuming after failures
for file in *.txt; do
```

### README Updates
When adding features:
1. Update relevant README.md
2. Add example usage
3. Document new requirements/dependencies
4. Update quick start if workflow changes

---

## Performance Considerations

### For Large Files
```bash
# Use streaming where possible
zcat large.gz | process | gzip > output.gz

# Not
zcat large.gz > temp
process temp > output
gzip output
```

### For Array Jobs
```bash
# Good: Process N samples in parallel
#SBATCH --array=1-100

# Check for existing output (idempotent)
if [[ -e "$output" ]]; then
    echo "Output exists, skipping"
    exit 0
fi
```

---

## Debugging Tips

### Enable Tracing
```bash
# At top of script
set -x  # Print each command before execution

# Or run script with
bash -x script.sh
```

### Check Variables
```bash
# Print all variables
set

# Print specific
echo "DEBUG: VAR=$VAR" >&2
```

### Validate Each Step
```bash
# Check critical outputs
if [[ ! -s "$output_file" ]]; then
    echo "Error: Output file missing or empty" >&2
    exit 1
fi
```

---

## Resources

- **ShellCheck Wiki**: https://www.shellcheck.net/wiki/
- **Bash Guide**: https://mywiki.wooledge.org/BashGuide
- **Google Shell Style Guide**: https://google.github.io/styleguide/shellguide.html
- **PEP 8**: https://pep8.org/
- **Python Type Hints**: https://docs.python.org/3/library/typing.html

---

## Contributing

When contributing to atavide_lite:

1. **Run shellcheck** on modified scripts
2. **Run black/flake8** on Python files
3. **Test on small dataset** before full run
4. **Update documentation** (this file, README, etc.)
5. **Clear commit messages** explaining what and why
6. **Consider backward compatibility** - don't break existing workflows

---

## CI/CD (Optional)

### GitHub Actions Example

If adding CI, keep it lightweight:

```yaml
# .github/workflows/lint.yml
name: Lint

on: [push, pull_request]

jobs:
  shellcheck:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Run shellcheck
        run: |
          sudo apt-get install -y shellcheck
          find . -name "*.sh" -o -name "*.slurm" | xargs shellcheck
        continue-on-error: true  # Don't fail builds, just warn

  python-lint:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: actions/setup-python@v4
        with:
          python-version: '3.10'
      - name: Install tools
        run: pip install black flake8
      - name: Check formatting
        run: black --check bin/*.py
        continue-on-error: true
```

**Keep CI non-blocking initially** - use as advisory, not gate.
