# atavide_lite `pawsey_minion` — agent runbook (Setonix / Pawsey)

Operational guidance for running the **long-read (ONT) single-end** atavide_lite pipeline on
Setonix. This file is deliberately generic: it records durable, reusable knowledge — decisions,
constraints and traps — not the log of any one project. Keep per-project state in that
project's own `.agent/CONTINUITY.md`, not here.

If you learn something durable while running this pipeline, add it here. If it is only true
for one dataset, it belongs in that project's continuity file.

---

## 1. Pick the right pipeline

| Your data | Use |
|---|---|
| ONT / PacBio, single-end | **`pawsey_minion/`** ← this directory |
| Illumina, paired `_R1`/`_R2` | `pawsey_shortread/` |

`pawsey_shortread/` requires paired files and **must not** be used with single-end long reads.
Never manufacture `_R1`/`_R2` filenames to satisfy it.

The source checkout is at **`$HOME/GitHubs/atavide_lite`** — not `$HOME/atavide_lite`. Several
scripts reference the checkout path; if you see the shorter form anywhere, it is a bug.

## 2. Preflight

```bash
cd <analysis directory>                       # must contain fastq/ and DEFINITIONS.sh
export SRC=$HOME/GitHubs/atavide_lite/pawsey_minion

find fastq -maxdepth 1 -type f -name '*.fastq.gz' -printf '%f\n' | sort > reads.txt
export NUM_READS=$(wc -l < reads.txt)

mkdir -p slurm_output/{host_slurm,megahit_slurm,mmseqs_slurm,vamb_slurm,fastplong_slurm}

pushd "$HOME/GitHubs/atavide_lite/bin" && make all && popd   # builds fastq2fasta, fastg2gfa, fasta_split
```

- `reads.txt` holds **one FASTQ filename per line** (basename only, no directory).
- **Create the `slurm_output/` directories before `sbatch`.** Slurm opens the job's stdout and
  stderr *before* the script runs, so a missing log directory fails the job with no useful
  message.
- `DEFINITIONS.sh` must set `SAMPLENAME`, `FILEEND`, `SOURCE`, `HOSTFILE`, `HOST`, `HOSTREMOVED`.
  Optionally `SPIKE_IN` and `SPIKE_IN_SEQUENCE` (see §5).
- Conda prefixes: `/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite` and
  `…/atavide_lite_vamb`. `pawsey_lib/check_atavide_lite_env.sh` validates the former.

### Verify databases exist before submitting anything

`/scratch` is subject to Pawsey's purge policy, so **databases silently disappear.** Check,
don't assume — a missing database is much cheaper to find now than three stages in:

```bash
for d in human UniRef50 UniRef100 NCBI/taxonomy; do
  p=/scratch/$PAWSEY_PROJECT/$USER/Databases/$d
  [ -e "$p" ] && echo "OK      $d ($(du -sh $p 2>/dev/null | cut -f1))" || echo "MISSING $d"
done
```

An empty-but-present directory counts as missing. Re-fetch with `download_human.slurm`,
`download_uniref50.slurm`, `download_uniref100.slurm`, `download_taxon_db.slurm`.

## 3. Slurm resources on Setonix — read this before changing any `--mem`

Setonix enforces **`MaxMemPerCPU`: 1840 MB on `work`, 7900 MB on `highmem`.**

A bare `--mem` larger than `cpus-per-task × MaxMemPerCPU` does **not** fail. Slurm silently
raises the allocated CPU count until it can satisfy the memory request. `--cpus-per-task=1
--mem=128G` therefore queues for **72 CPUs** — potentially hours of extra wait for a job that
uses one core for a couple of minutes.

These scripts request memory with **`--mem-per-cpu`** so the two figures can never disagree.
Keep it that way. If you change a memory request, change `--cpus-per-task` with it and confirm
what you actually got:

```bash
scontrol show job <jobid> | grep -E 'NumCPUs|ReqTRES'
```

`mmseqs_easy_taxonomy.slurm` and `download_uniref100.slurm` need more memory per core than
`work` provides, so they set `--partition=highmem` explicitly.

Derive tool thread counts from the allocation, never hardcode them:

```bash
THREADS=${SLURM_CPUS_PER_TASK:-8}
```

Hardcoded thread counts drift out of sync with `--cpus-per-task` and silently over-subscribe.
Watch for tools that auto-detect CPUs and see all 128 cores on the node rather than your
allocation — pass the thread count explicitly.

### Array throttling

Submit the whole array (`--array=1-$NUM_READS`) and let the scheduler decide concurrency. Add a
`%N` cap **only** for a concrete reason, and write the reason down. Legitimate reasons: an
external rate limit, contention for a resource Slurm doesn't manage, or a stage whose per-task
footprint is so large that unbounded concurrency would exhaust a shared filesystem. "Being
conservative" is not a reason.

## 4. Quality control — use `fastplong`, not `fastp`

**`fastplong.slurm` is the correct QC step for ONT data.** `fastp` is short-read software; its
upstream documentation points long-read users at fastplong.

This matters a great deal in practice. On a real ONT dataset, `fastp` with its default quality
filtering discarded **93% of reads** (1,014,189 of 1,089,169 failed `low_quality_reads`) — not
from adapter or length trimming, but because fastp judges a whole read by its global fraction
of low-quality bases. fastplong instead uses window-based `--break`, discarding low-quality
*regions* and keeping the rest of the read.

`fastp.slurm` remains only for backwards compatibility.

### The output-directory trap

`fastplong.slurm` writes to **`fastq_fastplong/`**, while `host_removal.slurm` defaults to
**`QC_DIR=fastq_fastp`**. So after a fastplong run you **must** pass:

```bash
--export=ALL,QC_DIR=fastq_fastplong
```

Without it, host removal either fails outright or — far worse — silently consumes stale `fastp`
output. This has caused real mixed-provenance incidents: in one run, two of sixteen samples
carried fastp-derived reads (4 MB) while the other fourteen were fastplong-derived (78 MB),
because a cancelled earlier attempt had left outputs behind. The size discrepancy was the only
visible clue.

**If you cancel a partially-completed stage, delete its derived outputs before re-running.**
`host_removal.slurm` refuses to overwrite existing spike-in outputs (exit 2) — that guard is
doing you a favour, so fix the provenance rather than working around it. Uniform provenance is
non-negotiable for any analysis someone will interpret.

`sankey_plot.slurm` and `read_fate.slurm` take `QC_DIR` and `QC_LABEL` so the QC transition is
labelled with the tool actually used.

### fastplong parameter notes

- **`-n` differs between the tools.** fastplong's short `-n` is an N *percentage*; fastp's is an
  absolute N-base count. Use the long form `--n_base_limit` to stay unambiguous.
- **Keep adapter trimming enabled.** Disabling it (`-A`) raises retention dramatically — 99.6%
  versus ~34% in one test — but that is retained adapter sequence, not recovered signal. Do not
  chase a retention target with `-A` unless library-level evidence shows the reads are
  adapter-free.
- Do not expect 80–90% retention from ONT metagenomes with adapter trimming on. Retention in the
  30–50% range is normal; forcing it higher means keeping adapter-dominated and very short reads.

## 5. Host and spike-in removal

`host_removal.slurm` maps with `minimap2 -x map-ont` against `$HOSTFILE` and keeps unmapped
reads. Peak RSS against GRCh38 + hs38d1 is around 11–12 GB, so 16 CPUs at 1800M/CPU is
comfortable.

If `DEFINITIONS.sh` sets `SPIKE_IN_SEQUENCE` (and optionally `SPIKE_IN`), the script then maps
host-removed reads to that reference, keeps only confident primary alignments (MAPQ ≥ 20), and
removes them. `$HOSTREMOVED` is therefore spike-free; the unfiltered intermediate is preserved
as `${HOSTREMOVED}_before_${SPIKE_IN}`.

**Check whether your run used a spike-in and declare it.** An undeclared spike will dominate
downstream tables — an MS2 process control accounted for 98.5% of reads in one negative-control
barcode. Record both host and spike-in counts; they are QC signal, not just waste.

## 6. Workflow

```bash
JOB=$(sbatch --parsable --array=1-$NUM_READS "$SRC/fastplong.slurm")

HOSTJOB=$(sbatch --parsable --array=1-$NUM_READS --dependency=afterok:$JOB \
    --export=ALL,QC_DIR=fastq_fastplong "$SRC/host_removal.slurm")

FAJOB=$(sbatch --parsable --dependency=afterok:$HOSTJOB "$SRC/fastq2fasta.slurm")

MMSEQSJOB=$(sbatch --parsable --array=1-$NUM_READS --dependency=afterok:$FAJOB \
    --export=ALL,MMSEQS_DB=UniRef50 "$SRC/mmseqs_easy_taxonomy.slurm")

MMTAXJOB=$(sbatch --parsable --dependency=afterok:$MMSEQSJOB \
    --export=ALL,MMSEQS_DB=UniRef50 "$SRC/mmseqs_summarise_taxonomy.slurm")
FATEJOB=$(sbatch --parsable --dependency=afterok:$MMSEQSJOB \
    --export=ALL,QC_DIR=fastq_fastplong,QC_LABEL=fastplong "$SRC/read_fate.slurm")
SSJOB=$(sbatch --parsable --array=1-$NUM_READS --dependency=afterok:$MMSEQSJOB \
    "$SRC/mmseqs_add_subsystems_taxonomy_fast.slurm")
COUNTSSJOB=$(sbatch --parsable --dependency=afterok:$SSJOB "$SRC/count_subsystems.slurm")
SANKEYJOB=$(sbatch --parsable --dependency=afterok:$COUNTSSJOB \
    --export=ALL,QC_DIR=fastq_fastplong,QC_LABEL=fastplong "$SRC/sankey_plot.slurm")

# Optional assembly / binning branch — independent of the read-based annotation above
MEGAHITJOB=$(sbatch --parsable --array=1-$NUM_READS --dependency=afterok:$HOSTJOB "$SRC/megahit.slurm")
VCJOB=$(sbatch --parsable --dependency=afterok:$MEGAHITJOB "$SRC/vamb_concat.slurm")
VMJOB=$(sbatch --parsable --array=1-$NUM_READS --dependency=afterok:$VCJOB "$SRC/vamb_minimap.slurm")
VAMBJOB=$(sbatch --parsable --dependency=afterok:$VMJOB "$SRC/vamb.slurm")
CHECKMJOB=$(sbatch --parsable --dependency=afterok:$VAMBJOB "$SRC/checkm.slurm" vamb/bins/ vamb/checkm)
```

### Choosing the MMseqs database

`mmseqs_easy_taxonomy.slurm` defaults to **UniRef100**; set `MMSEQS_DB` to override. Pass the
same value to `mmseqs_summarise_taxonomy.slurm`.

| | size | when |
|---|---|---|
| `UniRef50` | ~20 GB | first-pass and time-critical work — loads faster, more samples run concurrently |
| `UniRef100` | ~92 GB | final/deep analysis; more hits |

MMseqs against UniRef100 is by far the slowest stage — single tasks have run 19+ hours. Budget
walltime accordingly and prefer UniRef50 when speed matters.

### `--wait` versus dependencies

```
Does the agent need to inspect a result before deciding the next step?
    YES -> sbatch --wait
    NO  -> sbatch --parsable + afterok dependencies
```

Chain deterministic stages with dependencies rather than polling. Record job IDs.

## 7. Known traps

- **Completed arrays age out of the dependency table.** `--dependency=afterok:<jobid>` against an
  array that finished a while ago fails with `Batch job submission failed: Job dependency
  problem`. Bypass this only after explicitly verifying every task's exit code *and* the expected
  outputs, then submit a fresh chain with its own internal dependencies.
- **Incomplete trailing FASTQ records** (a truncated final record) break record counting.
  `sankey_plot.py` logs and skips incomplete trailing records but still errors on malformed
  complete ones. Investigate rather than assuming truncation is benign.
- **`vamb.slurm` requests a GPU** and hardcodes a `-gpu` account. Don't add an account override
  unless cluster policy requires it. GPU jobs need the `${PAWSEY_PROJECT}-gpu` account.
- **Never run compute on a login node.** Databases here are tens to hundreds of GB; loading one
  on a shared login node degrades it for every other user and will likely be OOM-killed anyway.
- **Don't re-run into populated output directories.** Use a fresh analysis directory (or a new
  output namespace) so a comparison run cannot clobber or half-overwrite prior results.
- **`mmseqs_easy_taxonomy.slurm` skips a sample whose output directory already exists**
  (`"$OUTPUT exists. Nothing to do"`, then `exit 0`). This is a **silent** failure mode after a
  cancelled or failed run: the previous attempt leaves behind an empty `mmseqs/<sample>/`
  directory, the rerun sees it, skips the sample, and **reports success**. Any `afterok`
  dependency is satisfied and the pipeline proceeds with that sample missing entirely.

  Observed in practice: a cancelled task left an empty output directory; on resubmission the
  task "COMPLETED 0:0" in seconds and one of sixteen samples had no taxonomy at all.

  **After cancelling any stage, delete the output directories of the affected tasks before
  resubmitting**, and sanity-check that each task actually produced output rather than trusting
  its exit code:

  ```bash
  # a task that really ran leaves files; a skipped one leaves an empty directory
  find mmseqs -mindepth 1 -maxdepth 1 -type d -empty
  ```

  More generally on this pipeline: **exit code 0 is not proof a stage did the work.** Always
  confirm the expected outputs exist and are non-empty.

## 8. Validation before submitting

```bash
bash -n <script>.slurm                 # syntax
python3 -c "import ast,sys; ast.parse(open(sys.argv[1]).read())" <script>.py
git diff --check                       # whitespace, if you edited the repo
```

Then confirm the allocation matches intent with `scontrol show job <jobid>`.

## 9. Reporting

Report materially: job IDs and arrays submitted, whether they completed, meaningful failures and
what was changed in response, resource observations, and output locations. Distinguish a
verified success from an assumed one — never report a stage as complete without checking exit
codes *and* that the expected outputs exist.
