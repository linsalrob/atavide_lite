# Known-Good Software Versions

This document lists tested versions of tools used in the atavide_lite pipeline. While the pipeline may work with other versions, these represent tested baselines for reproducibility.

---

## Core Pipeline Tools

| Tool | Tested Version | Notes |
|------|----------------|-------|
| **fastp** | v0.23.2+ | Quality trimming and adapter removal |
| **minimap2** | v2.29+ | Read alignment for host removal |
| **samtools** | v1.20+ | SAM/BAM processing |
| **mmseqs2** | latest | Protein sequence search and taxonomy assignment |
| **megahit** | v1.2.9+ | Metagenomic assembly |
| **vamb** | v4.0.0+ | Binning tool; **see compatibility notes** below |
| **checkm-genome** | v1.2.0+ | MAG quality assessment |

---

## Supporting Tools

| Tool | Tested Version | Notes |
|------|----------------|-------|
| **Python** | 3.9 - 3.11 | Required for helper scripts |
| **pigz** | latest | Parallel gzip compression |
| **parallel** | latest | GNU parallel for batch processing |
| **snakemake** | v7.0+ | Used for some sub-workflows |
| **pytaxonkit** | latest | Taxonomic utilities |
| **taxonkit** | latest | Taxonomy data manipulation |
| **rclone** | latest | Remote file transfer (optional) |
| **rsync** | latest | File synchronization |
| **sra-tools** | v3.0.0+ | SRA data download (optional) |

---

## Database Versions and Snapshots

Databases change over time. Recording the version/date of databases used ensures reproducibility.

### UniRef (Protein Database)

- **Database**: UniRef50 or UniRef90
- **How to record version**: 
  ```bash
  # Check database creation date
  ls -lh /path/to/mmseqs/db/*.dbtype
  # Save download date
  echo "Downloaded: $(date +%Y-%m-%d)" > UniRef50_version.txt
  ```
- **Recommended**: Download every 6-12 months for updated annotations
- **Size**: UniRef50 ≈ 20-30 GB (compressed), UniRef90 ≈ 50-80 GB

### BV-BRC (Bacterial and Viral Bioinformatics Resource Center)

- **Mapping files**: Subsystems annotations
- **How to record**: Save download date and source URL
  ```bash
  # Example
  wget https://www.bv-brc.org/api/genome_feature/?eq(annotation,PATRIC) \
    -O bvbrc_snapshot_$(date +%Y-%m-%d).json
  ```
- **Update frequency**: Quarterly updates recommended

### Host Reference Genomes

| Host | Reference Assembly | Source |
|------|-------------------|--------|
| Human | GRCh38 (GCA_000001405.15) | NCBI/GenBank |
| Mouse | GRCm39 (GCA_000001635.9) | NCBI/GenBank |
| Other | Species-specific | Check NCBI, Ensembl, or specific databases |

**Recording version**:
```bash
# Example for human reference
grep '^>' GCA_000001405.15*.fna | head -1
# Save assembly name and download date
```

---

## Python Package Versions

For reproducibility of helper scripts (`bin/` directory):

```txt
# From requirements.txt (if available) or conda environment
biopython>=1.79
numpy>=1.21
pandas>=1.3.0
matplotlib>=3.4.0
seaborn>=0.11.0
plotly>=5.0.0  # For Sankey plots
pysam>=0.19.0
```

To export your exact environment:
```bash
conda env export > atavide_lite_environment.yml
# or
pip freeze > requirements_frozen.txt
```

---

## Cluster/Scheduler Software

| System | Scheduler | Tested Version |
|--------|-----------|----------------|
| Pawsey Setonix | Slurm | v22+ |
| Flinders Deepthought | Slurm | v21+ |
| NCI Gadi | PBS Pro | v19+ |

---

## VAMB Compatibility Notes

**CRITICAL**: VAMB has had breaking changes between versions.

### Known Issues:
- **VAMB v4.0.0+**: Changed argument structure
  - Removed or renamed some parameters (e.g., `minsize` handling)
  - API changes in `vamb.vambtools` module
- **VAMB v3.x**: Different cluster output format

### Recommended Approach:
1. **Pin VAMB version** in your conda environment:
   ```bash
   conda install vamb=4.1.0  # or specific tested version
   ```
2. **Document your version**:
   ```bash
   python -c "import vamb; print(vamb.__version__)" > vamb_version.txt
   ```
3. See `docs/compat.md` for version-specific handling in scripts

---

## How to Use This Document

### For a New Analysis:
1. Record versions of all tools at the start:
   ```bash
   fastp --version
   minimap2 --version
   samtools --version
   mmseqs --version
   megahit --version
   python -c "import vamb; print(vamb.__version__)"
   ```
2. Save database download dates and sources
3. Document any deviations from these tested versions

### For Reproducing Results:
1. Install tools matching these versions (use conda/mamba for exact versions)
2. Download databases from the same time period or use archived versions
3. Check for known incompatibilities if versions differ

### For Reporting Issues:
Include tool versions in bug reports to help diagnose environment-specific problems.

---

## Installation Commands

### Create Conda Environment with Pinned Versions:
```bash
# Using the provided YAML file
mamba env create -f atavide_lite.yaml

# Or manual installation with specific versions
mamba create -n atavide_lite python=3.10
mamba activate atavide_lite
mamba install fastp=0.23.2 minimap2=2.29 samtools=1.20 \
    mmseqs2 megahit=1.2.9 vamb=4.1.0 checkm-genome \
    snakemake=7.32.4 pytaxonkit taxonkit pigz parallel
```

### Verify Installation:
```bash
# Run this script to check versions
for tool in fastp minimap2 samtools mmseqs megahit checkm; do
    echo -n "$tool: "
    $tool --version 2>&1 | head -1
done

python -c "import vamb; print(f'vamb: {vamb.__version__}')"
```

---

## Notes

- **Not strict requirements**: The pipeline may work with other versions, but these are tested
- **Update this document**: When you test new versions, please update this file
- **Compatibility**: Newer versions usually work, but breaking changes can occur (especially VAMB)
- **HPC modules**: If using system modules instead of conda, record module versions:
  ```bash
  module list > module_versions_$(date +%Y-%m-%d).txt
  ```

---

## Reproducibility Best Practices

1. **Document everything**: Tool versions, database dates, parameters
2. **Use environment files**: Export conda environments for exact replication
3. **Archive databases**: Keep copies of database snapshots if possible
4. **Version control configs**: Track `DEFINITIONS.sh` and config files in git
5. **Record metadata**: Sample processing dates, batch information, any failures/reruns
