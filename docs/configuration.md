# Configuration and input preparation

## Analysis directory layout

Run the workflow from a dedicated analysis directory, not from the repository checkout. A typical paired-read analysis begins as:

```text
analysis/
├── DEFINITIONS.sh
├── R1_reads.txt
├── fastq/
│   ├── sample_a_R1.fastq.gz
│   └── sample_a_R2.fastq.gz
└── slurm_output/
```

Single/long-read analyses use `reads.txt`. Pipeline output directories are created alongside these control files. Keeping the repository and analysis separate prevents accidental commits of data and allows one code checkout to serve multiple analyses.

## `DEFINITIONS.sh`

Every scheduler job runs in the analysis directory and sources `DEFINITIONS.sh`. Copy an example from the selected profile and review every value.

`SAMPLENAME`
: A short dataset identifier used in output names. Avoid whitespace and shell metacharacters.

`FILEEND`
: The exact suffix removed to derive sample names. In a paired profile it normally includes the R1 marker, such as `_R1.fastq.gz`; in a single-read profile it may be `.fastq.gz`.

`FAFILEEND`
: Used by profiles that derive names from converted FASTA files. It generally mirrors `FILEEND` with the FASTQ extension replaced by FASTA. Not every current profile defines it.

`SOURCE`
: Directory containing the starting FASTQ files, commonly `fastq`.

`HOSTFILE`
: FASTA or compressed FASTA containing the reference sequence(s) to separate from downstream reads. It may be omitted when host removal is disabled, depending on the chosen definitions example.

`HOST`
: Output directory for reads mapped to the host reference.

`HOSTREMOVED`
: Output directory for reads retained for downstream metagenomics analysis. When host removal is skipped, profiles commonly point this at the fastp output directory.

Cluster profiles also rely on system variables such as `PAWSEY_PROJECT`, `USER`, `BGFS`, scheduler job/array IDs, or `PBS_JOBFS`. Confirm these are defined inside batch jobs rather than only in an interactive login shell.

## Paired-read naming

Paired scripts commonly derive R2 from R1 by replacing `_R1` with `_R2`. Each filename must therefore contain the expected marker exactly once. Rename `_1`/`_2` downloads if required:

```bash
for file in fastq/*_1.fastq.gz; do
  mv "$file" "${file/_1.fastq.gz/_R1.fastq.gz}"
done
for file in fastq/*_2.fastq.gz; do
  mv "$file" "${file/_2.fastq.gz/_R2.fastq.gz}"
done
```

Preview rename operations on valuable data and follow local data-management practice before modifying originals.

Create the sample list deterministically:

```bash
find fastq -type f -name '*_R1*' -printf '%f\n' | sort > R1_reads.txt
```

Check counts and mates:

```bash
wc -l R1_reads.txt
while IFS= read -r r1; do
  r2=${r1/_R1/_R2}
  test -f "fastq/$r2" || printf 'Missing mate: %s\n' "$r2"
done < R1_reads.txt
```

## Single/long-read naming

Create one entry per file:

```bash
find fastq -type f -name '*.fastq.gz' -printf '%f\n' | sort > reads.txt
```

Choose `FILEEND` so removing it yields a unique, meaningful sample name. Do not include reports, checksums, or intermediate FASTQ files in `reads.txt`.

## Host-removal choices

Host filtering is a scientific and governance decision. Select an appropriate reference assembly, record its accession and version, and understand whether removal is required before data leave controlled storage. A multi-FASTA can represent a non-human host, vector, symbiont, or other sequence collection.

The pipeline separates matching and non-matching reads, but downstream stages currently use the non-host set. Review mapping statistics and a small sample of reads before discarding or restricting access to either output.

## Database provenance

Record the following with each analysis:

- database name and release, such as UniRef100;
- download or build date and source URL;
- taxonomy dump date;
- Subsystems mapping/database version;
- host reference accessions and checksums; and
- commands used to create local MMseqs2 indices.

Changing a database changes the analysis even when the scripts and input reads are identical.

## Resource settings

Resource headers are starting points. Before submission, review `--time`, `--cpus-per-task`, `--mem`, nodes, partition/queue, QoS, account/project, array limits, output paths, and GPU directives in every selected script. Ensure the application's own thread option matches the scheduler CPU allocation; requesting one CPU while starting many threads creates contention and misleading accounting.
