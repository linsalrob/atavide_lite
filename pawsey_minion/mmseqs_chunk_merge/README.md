# Chunked MMseqs taxonomy on Pawsey

This opt-in workflow replaces one large `mmseqs easy-taxonomy` job per sample
with base-balanced FASTA chunks, an unthrottled Slurm search array, and one
validated merge job per sample. It leaves `mmseqs_easy_taxonomy.slurm`
unchanged.

The merge recreates all four canonical `easy-taxonomy` outputs under
`mmseqs/<sample>/`:

```text
<sample>_lca.tsv.gz
<sample>_report.gz
<sample>_tophit_aln.gz
<sample>_tophit_report.gz
```

Chunk FASTAs are complete at record boundaries and balanced by sequence bases.
Every source identifier must be unique and occur in exactly one emitted chunk.
Per-query LCA and alignment files are concatenated because chunk query sets are
disjoint. Aggregate taxonomy reports are regenerated from merged MMseqs result
databases; they are never concatenated. Top-hit reports are regenerated from
retained alignment databases and checked against additive chunk counts.

## Run

Run from an analysis directory containing `DEFINITIONS.sh`, `reads.txt`, and
the `fasta/` outputs from `fastq2fasta.slurm`:

```bash
export SRC="$HOME/GitHubs/atavide_lite/pawsey_minion/mmseqs_chunk_merge"
bash "$SRC/submit.sh"
```

By default, the splitter aims for approximately one billion uncompressed bases
per chunk. Override this and other paths without editing the scripts:

```bash
TARGET_BASES=750000000 \
WORK_ROOT=mmseqs_chunking \
OUTPUT_ROOT=mmseqs \
bash "$SRC/submit.sh"
```

`submit.sh` prints the split and preparation job IDs. The preparation job
submits the complete chunk array without a `%N` throttle, gives each sample
merge an `afterok` dependency on exactly its own chunks, and records all IDs in
`$WORK_ROOT/job_ids.tsv`. A pipeline is finished only when all merge jobs have
completed successfully and each sample has `merged/MERGE.DONE` plus the four
gzip-valid canonical outputs.

Defaults are 48 CPUs/~84 GiB/24 hours per chunk and 72 CPUs/~127 GiB/24 hours per
sample merge. The CPU counts look large for the work they do: on Setonix, `MaxMemPerCPU`
is 1840 MB on `work`, so CPUs are the only way to request memory and these counts are simply
what that much memory costs. Requesting it via `--mem-per-cpu` means Slurm allocates exactly
this rather than silently inflating a bare `--mem`. Resource directives can be overridden with ordinary `sbatch`
options if a dataset needs different allocations. The MMseqs executable and
UniRef100 database can be overridden with `MMSEQS_BIN` and `TARGET_DB`.
