# Quick start

This walkthrough uses paired-end data and Slurm-like commands to explain the shared workflow. Use the exact README and scheduler syntax in the selected profile. Start with one or two small, non-sensitive samples; do not begin with an entire production dataset.

## 1. Install and choose a profile

Complete [Installation](installation.md), then point `ATAVIDE_PROFILE` at the chosen directory. For Pawsey paired reads:

```bash
export ATAVIDE_PROFILE="$HOME/GitHubs/atavide_lite/pawsey_shortread"
```

Avoid naming this variable `SOURCE`: `SOURCE` is a dataset setting inside `DEFINITIONS.sh`.

## 2. Create an analysis directory

Keep input data and analysis outputs outside the cloned repository:

```bash
mkdir -p /path/to/analysis/fastq
cd /path/to/analysis
```

Place or link paired files into `fastq/` so every sample has one R1 and one R2 file, for example:

```text
fastq/
├── sample_a_R1.fastq.gz
├── sample_a_R2.fastq.gz
├── sample_b_R1.fastq.gz
└── sample_b_R2.fastq.gz
```

## 3. Configure the dataset

Copy the closest definitions example:

```bash
cp "$ATAVIDE_PROFILE/DEFINITIONS_human.sh" DEFINITIONS.sh
```

Edit it for the dataset. At minimum check `SAMPLENAME`, `FILEEND`, `SOURCE`, `HOSTFILE`, `HOST`, and `HOSTREMOVED`. For analysis without host removal, use a no-host example when supplied or set `HOSTREMOVED` to the quality-controlled read directory as documented by the profile.

```{warning}
`DEFINITIONS.sh` is sourced as shell code by jobs. Do not insert untrusted commands, spaces around assignments, or secrets. Quote paths in new adaptations and avoid whitespace in sample identifiers.
```

## 4. Build the sample list

For paired-end reads:

```bash
find fastq -name '*_R1*' -printf '%f\n' | sort > R1_reads.txt
export NUM_R1_READS="$(wc -l < R1_reads.txt)"
```

Check the list and verify every R1 has a corresponding R2. For single/long reads, create `reads.txt` instead and export `NUM_READS`:

```bash
find fastq -type f -name '*.fastq.gz' -printf '%f\n' | sort > reads.txt
export NUM_READS="$(wc -l < reads.txt)"
```

These counts are shell-session variables. Recreate them after logging in again or store a safe setup command with the analysis notes.

## 5. Create log directories

Profiles write scheduler logs to stage-specific paths. Create those used by the selected scripts before submission:

```bash
mkdir -p \
  slurm_output/fastp_slurm \
  slurm_output/host_slurm \
  slurm_output/megahit_slurm \
  slurm_output/mmseqs_slurm \
  slurm_output/vamb_slurm
```

## 6. Prepare databases

Confirm the host reference, UniRef/MMseqs2 database, taxonomy data, Subsystems mapping, and optional tool databases exist at the paths expected by the profile. On Pawsey, profile download scripts can be submitted separately. Wait for required databases before submitting dependent annotation jobs.

## 7. Submit quality control

For paired Slurm profiles:

```bash
QC_JOB=$(sbatch --parsable \
  --array="1-${NUM_R1_READS}:1" \
  "$ATAVIDE_PROFILE/fastp.slurm")
echo "Quality-control job: $QC_JOB"
```

Inspect the first completed `.out` and `.err` files and fastp reports before continuing.

## 8. Submit optional host removal

```bash
HOST_JOB=$(sbatch --parsable \
  --array="1-${NUM_R1_READS}:1" \
  --dependency="afterok:${QC_JOB}" \
  "$ATAVIDE_PROFILE/host_removal.slurm")
```

If host removal is disabled, identify the cleaned-read job and directory that downstream scripts should use; do not blindly submit a nonexistent dependency.

## 9. Start the two analysis branches

Assembly can begin from cleaned reads without waiting for read-based taxonomy:

```bash
ASSEMBLY_JOB=$(sbatch --parsable \
  --dependency="afterok:${HOST_JOB}" \
  "$ATAVIDE_PROFILE/megahit_allreads.slurm")
```

Prepare FASTA for MMseqs2:

```bash
FASTA_JOB=$(sbatch --parsable \
  --dependency="afterok:${HOST_JOB}" \
  "$ATAVIDE_PROFILE/fastq2fasta.slurm")
```

Then submit the taxonomy script available in the selected profile. Pawsey paired profiles currently provide UniRef50 and UniRef100 variants; choose one deliberately and keep the resulting database identity with the analysis metadata.

## 10. Monitor before scaling up

Use the scheduler's status and accounting commands to check exit status, elapsed time, maximum memory, and CPU efficiency. Inspect outputs after every new stage on the small dataset. Only then increase array size or process the full dataset.

Continue with [Pipeline stages](pipeline.md) for downstream taxonomy, function, assembly, VAMB, and CheckM dependencies.
