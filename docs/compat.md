# Compatibility Documentation

This document tracks known compatibility issues and how we handle them in the atavide_lite pipeline.

---

## VAMB Version Compatibility

### The Problem

VAMB has undergone significant API changes across versions, causing scripts that work with one version to fail with another. This is a common issue in rapidly developing tools.

### Known Breaking Changes

#### VAMB v4.0.0+ Changes:
1. **Parameter changes in bin creation**:
   - `minsize` parameter behavior changed
   - Some API functions were renamed or restructured
   - Module structure changed (`vamb.vambtools` interface)

2. **Cluster file format**:
   - Earlier versions: Different TSV structure
   - v4.0.0+: Standardized `clustername<TAB>contigname` format

3. **Python API**:
   - `vamb.vambtools.read_clusters()` behavior changed
   - `vamb.vambtools.write_bins()` signature updated

### Our Solution

We maintain **two VAMB helper scripts** with different compatibility approaches:

#### 1. `bin/vamb_create_fasta_clusters.py` (Canonical)
- **Status**: ✅ **Primary/Canonical script**
- **Dependencies**: Minimal - only standard library + atavide_lib
- **Compatibility**: Works with any VAMB version (doesn't import vamb)
- **Usage**: Reads cluster TSV and creates FASTA files manually
- **Pros**: 
  - Version-independent
  - No VAMB import required
  - Easier to debug
  - Works with both v3.x and v4.x cluster formats
- **Cons**: 
  - Slightly slower for very large datasets
  - Reimplements some VAMB functionality

#### 2. `bin/vamb_create_fasta.py` (Deprecated)
- **Status**: ⚠️ **Deprecated - kept for compatibility**
- **Dependencies**: Requires vamb module
- **Compatibility**: Designed for VAMB v4.0.0+
- **Issue**: Breaks with other VAMB versions due to API changes
- **Recommendation**: Use `vamb_create_fasta_clusters.py` instead

### Version Detection Pattern

If you need to write VAMB-dependent code, use this pattern:

```python
#!/usr/bin/env python3
import sys

# Attempt to detect VAMB version
try:
    import vamb
    try:
        VAMB_VERSION = vamb.__version__
    except AttributeError:
        # Older versions may not have __version__
        VAMB_VERSION = "unknown"
except ImportError:
    print("Error: VAMB not installed", file=sys.stderr)
    print("Install with: conda install vamb", file=sys.stderr)
    sys.exit(1)

# Check version compatibility
if VAMB_VERSION != "unknown":
    major_version = int(VAMB_VERSION.split('.')[0])
    if major_version < 4:
        print(f"Warning: VAMB {VAMB_VERSION} detected. Tested with VAMB 4.x", 
              file=sys.stderr)
        print("Consider upgrading: conda install vamb>=4.0.0", file=sys.stderr)
else:
    print("Warning: Could not determine VAMB version", file=sys.stderr)
    print("This script is tested with VAMB 4.0.0+", file=sys.stderr)

# Version-specific logic
if VAMB_VERSION.startswith('4.'):
    # Use v4 API
    from vamb.vambtools import read_clusters, write_bins
else:
    # Use v3 API or fail gracefully
    print("Error: Unsupported VAMB version", file=sys.stderr)
    sys.exit(1)
```

### Recommended Workflow

1. **Use the canonical script**: `vamb_create_fasta_clusters.py`
   - This avoids VAMB version issues entirely
   - Works consistently across environments

2. **Pin VAMB version** in conda environment:
   ```bash
   conda install vamb=4.1.0
   ```

3. **Document VAMB version** in your analysis:
   ```bash
   python -c "import vamb; print(vamb.__version__)" > vamb_version.txt
   ```

---

## MMseqs2 Database Compatibility

### The Problem

MMseqs2 databases have version-specific binary formats. Databases created with one version may not be readable by another.

### Solution

- **Always rebuild databases** when upgrading MMseqs2 major versions
- **Document database MMseqs2 version**:
  ```bash
  mmseqs version > mmseqs_db_build_version.txt
  ```
- **Keep databases separate by version**: 
  ```
  /databases/
    mmseqs2_v14/
      UniRef50/
    mmseqs2_v15/
      UniRef50/
  ```

---

## Python Version Compatibility

### Minimum Version: Python 3.9

#### Issues with older Python:
- Type hints (`dict[str, int]`) require Python 3.9+
- Some dependencies (numpy, pandas) dropping 3.8 support

#### Issues with very new Python:
- Some HPC systems lag in Python versions
- Conda packages may not be available immediately for Python 3.12+

### Recommendation
- **Target**: Python 3.9-3.11
- **Test new scripts** with both 3.9 and 3.11

---

## Conda vs System Modules

### The Problem

HPC systems often provide modules (via `module load`) that may conflict with conda environments.

### Best Practice

1. **Choose one approach per session**:
   ```bash
   # Option 1: Conda only
   module purge
   conda activate atavide_lite
   
   # Option 2: System modules only
   module load fastp/0.23.2
   module load minimap2/2.29
   # etc.
   ```

2. **Document which approach**:
   - Add to logs: `module list` or `conda list`
   - Include in methods section of papers

3. **Avoid mixing**:
   - Don't `module load samtools` then use conda's minimap2
   - PATH conflicts can cause subtle errors

---

## Slurm vs PBS

### Key Differences

| Feature | Slurm | PBS |
|---------|-------|-----|
| Array variable | `$SLURM_ARRAY_TASK_ID` | `$PBS_ARRAYID` |
| Job ID | `$SLURM_JOB_ID` | `$PBS_JOBID` |
| Temp storage | `$TMPDIR` | `$PBS_JOBFS` (NCI) |
| Array jobs | Supported | Not on NCI Gadi |

### Our Approach

- Maintain separate scripts for Slurm and PBS systems
- Share common logic via `lib/common.sh` (sourced by all)
- Document system-specific variables in script comments

---

## Filesystem Differences

### Fast Local Storage

Different HPC systems have different fast local storage:

| System | Fast Storage | Usage |
|--------|-------------|-------|
| Pawsey Setonix | `/scratch` | Shared across nodes |
| Flinders Deepthought | `$BGFS` | Node-local (not shared) |
| NCI Gadi | `$PBS_JOBFS` | Job-specific temporary |

### Pattern for Portable Scripts

```bash
# Detect system and set fast storage
if [[ -n "$PBS_JOBFS" ]]; then
    # NCI Gadi
    FAST_STORAGE="$PBS_JOBFS"
elif [[ -n "$BGFS" ]]; then
    # Flinders Deepthought
    FAST_STORAGE="$BGFS"
elif [[ -d "/scratch" ]]; then
    # Pawsey or similar
    FAST_STORAGE="/scratch/${PROJECT}/${USER}/tmp"
else
    # Fallback
    FAST_STORAGE="${TMPDIR:-/tmp}"
fi
```

---

## Recording Compatibility Issues (Template)

When you encounter a new compatibility issue:

1. **Document it here** with this structure:

```markdown
## [Tool Name] [Issue Description]

### The Problem
[Describe what breaks and why]

### Affected Versions
- Works: [version]
- Breaks: [version]

### Our Solution
[How we handle it in scripts]

### Workarounds
[Alternative approaches if needed]

### Test Case
[How to reproduce/test the issue]
```

2. **Update scripts** with version checks or fallback logic

3. **Note in commit message**: "Fix compatibility: [brief description]"

---

## Future Compatibility Considerations

### When adding new tools:
1. Document tested version in `docs/known_good_versions.md`
2. Add version detection if API is unstable
3. Fail gracefully with helpful error messages
4. Consider version-independent approaches (like our VAMB solution)

### When updating tools:
1. Test on small dataset first
2. Check for API changes in release notes
3. Update compatibility notes here
4. Update version in documentation

---

## Getting Help

If you encounter a compatibility issue:

1. **Check versions**:
   ```bash
   tool --version
   python -c "import module; print(module.__version__)"
   ```

2. **Check this document** for known issues

3. **Search GitHub issues** for the tool

4. **Test with known-good versions** from `docs/known_good_versions.md`

5. **Open an issue** in atavide_lite repo with:
   - Tool versions
   - Error messages
   - System details (HPC system, scheduler, etc.)
