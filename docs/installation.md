# Installation

## Prerequisites

Use a Linux HPC environment with Bash, Git, a C compiler and `make`, and either Slurm or PBS for the supplied cluster profiles. Conda-compatible package management is recommended; the examples use Mamba because dependency resolution is usually faster.

Before installing, confirm that the intended compute nodes can access:

- the cloned repository;
- the working and database directories;
- the selected temporary filesystem; and
- the Conda environment or modules used by batch jobs.

Do not perform computational analysis on a login node. Compilation and environment creation should also follow local HPC policy.

## Clone the repository

The profile scripts commonly expect the repository under `$HOME/GitHubs`, so that layout requires the fewest changes:

```bash
mkdir -p "$HOME/GitHubs"
cd "$HOME/GitHubs"
git clone https://github.com/linsalrob/atavide_lite.git
cd atavide_lite
```

To record the exact code used in an analysis:

```bash
git rev-parse HEAD
```

Store that commit ID with the analysis metadata. Pulling later changes during an active analysis can make stages inconsistent, so update deliberately.

## Build the bundled utilities

Compile the C utilities in `bin/`:

```bash
cd "$HOME/GitHubs/atavide_lite/bin"
make all
```

This builds at least `fastq2fasta` and `fastg2gfa`. Confirm the build completed without compiler errors before submitting dependent jobs.

## Create software environments

The main environment contains tools used across the read-processing and annotation stages:

```bash
mamba env create \
  --file "$HOME/GitHubs/atavide_lite/atavide_lite.yaml"
```

VAMB has a separate environment definition:

```bash
mamba env create \
  --file "$HOME/GitHubs/atavide_lite/atavide_lite_vamb.yaml"
```

Important dependencies include fastp, minimap2, samtools, MMseqs2, MEGAHIT, VAMB, CheckM, TaxonKit/PyTaxonKit, Snakemake, GNU Parallel, pigz, rsync, and rclone. The environment files are the source of truth for packaged dependencies.

```{note}
The Pawsey scripts intentionally place environments under `/scratch/$PAWSEY_PROJECT/$USER` and provide `pawsey_lib/check_atavide_lite_env.sh` because scratch content is temporary. Follow the Pawsey profile README instead of assuming a named environment in `$HOME`.
```

## Databases and references

Pipeline stages may require:

- a host reference FASTA for optional host removal;
- an MMseqs2 UniRef database, with current profiles favouring UniRef100 where provided;
- NCBI taxonomy dumps for TaxonKit;
- the atavide lite BV-BRC Subsystems SQLite mapping; and
- CheckM or PHROG data when those optional stages are used.

Pawsey profiles include scheduler scripts for several database downloads. On other systems, place databases on storage visible to compute nodes and update the profile paths. Databases are large, change over time, and affect reproducibility; record their source, release/version, retrieval date, and local path.

## Verify before using real data

Run these checks from the repository:

```bash
test -x bin/fastq2fasta
test -x bin/fastg2gfa
bash -n pawsey_shortread/fastp.slurm
```

Then submit a minimal scheduler job that prints the hostname and tool versions. Only after storage, environment activation, modules, accounting, and logging work should you run a small non-sensitive sample through the pipeline.

## Updating

Check the current branch and local changes before updating:

```bash
git status --short --branch
git pull --ff-only
```

Rebuild `bin/` and review environment changes after an update. Do not overwrite locally adapted cluster scripts without first saving them on their own branch.
