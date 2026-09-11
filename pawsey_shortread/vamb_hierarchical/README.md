# Hierarchical VAMB clustering

This directory contains a scalable, stand-alone continuation of VAMB from its
feature checkpoint. It is intended for datasets where conventional global
VAMB clustering cannot finish within the scheduler walltime or where writing a
very large number of bins makes one monolithic job difficult to recover.

The workflow starts from a completed VAMB feature directory containing:

- `composition.npz`
- `abundance.npz`
- `latent.npz`

It does not modify VAMB itself or the original feature files. Its primary
product is a conventional `clustername<TAB>contigname` assignment table, from
which normal FASTA-bin and quality-control tools can continue.

## Why a hierarchy?

VAMB's `ClusterGenerator` is a destructive, greedy medoid algorithm. Every
accepted bin changes the set of contigs available to the next decision, so its
global clustering loop is inherently serial. On very large latent matrices,
moving that loop to a GPU is not enough to guarantee completion inside a fixed
walltime. Splitting each individual distance calculation across GPUs can also
cost more in device transfers and coordination than it saves.

This workflow instead changes the unit of parallelism:

```text
VAMB feature checkpoint
        |
        v
GPU partition of latent space into S spatially coherent shards
        |
        v
S independent local VAMB jobs (Slurm GPU array)
        |
        v
CPU collection of local-bin representatives
        |
        v
global VAMB reconciliation of multi-contig representatives (GPU)
        |
        v
expand to contig assignments -> FASTA bins -> CheckM
        |
        v
optional stability-gated split proposals for suspicious bins
```

The number of shards is configurable; it is not tied to the number of GPUs.
Slurm schedules one GPU per local array task and runs as many tasks
concurrently as the account and cluster can support. More shards make the
local problems smaller and easier to recover, but increase scheduling,
filesystem, and boundary-reconciliation overhead. A value such as 256 is a
reasonable starting point for a very large dataset, not a universal default.

This method is an approximation to monolithic VAMB clustering. The hierarchy
makes the approximation explicit, retains provenance through every stage, and
adds a global reconciliation step to repair likely shard-boundary splits.

## The clustering stages

### 1. Partition the VAMB embedding

`hierarchical_vamb_clustering_rocm.py partition` applies spherical k-means to
the normalized VAMB latent coordinates on one GPU. The deterministic partition
creates balanced-enough, spatially coherent shard membership files. Shards
describe work units, not biological bins.

### 2. Cluster shards independently

A Slurm array runs the normal VAMB clustering algorithm within every shard.
Each task requests one GPU; no fixed GPU count is assumed. The array is left
unthrottled by default so Slurm and the account limits determine safe
concurrency.

Local tasks are restartable. A shard is complete only when both
`local_clusters.tsv` and `representatives.npz` exist. Outputs are first written
to temporary files and atomically renamed, so a timeout or filesystem problem
cannot leave a one-sided partial result that is mistaken for success. Failed
array elements can therefore be resubmitted independently.

### 3. Represent local bins

The CPU collection stage replaces each local bin with a compact record: its
latent centroid, total assembled bases, contig count, and mapping back to the
local members. This greatly reduces the number of objects presented to the
global stage.

Local singleton bins are kept explicitly. Only multi-contig representatives
need global clustering; preserving singletons separately matches VAMB's loner
semantics and prevents contigs from disappearing during reconciliation.

### 4. Reconcile across shard boundaries

The CPU prepare step builds the candidate representative matrix. A GPU job then
runs VAMB over the multi-contig representatives, allowing compatible local bins
from different shards to merge. The CPU expansion step translates those
representative clusters back to original contigs and combines them with the
retained singleton representatives.

The canonical result is:

```text
vamb_crass/SAMPLE/hierarchical_vamb/hierarchical_unsplit.tsv
```

FASTA export creates:

```text
vamb_crass/SAMPLE/hierarchical_vamb/hierarchical_unsplit_bins/
```

The exporter uses the same minimum-bin-size convention as a normal VAMB
workflow. The assignment table remains the authoritative complete mapping;
the FASTA directory may omit bins below the export threshold.

### 5. CheckM and optional refinement proposals

CheckM runs on CPU and does not require a GPU. The supplied launchers put its
large per-bin intermediates on node-local temporary storage and retain the
tabular `checkm/qa.tsv` result. This avoids exhausting a shared filesystem's
file-count quota.

For very large bin collections, the batch workflow prepares an absolute-path
manifest, runs independent CheckM array batches, and merges their QA tables.
Batching changes scheduling and the failure domain, not the biological
analysis. Individual failed batches can be retried without repeating all bins.

CheckM can identify bins that may have been over-merged, but contamination
alone is not sufficient evidence to split a bin. The current gated evaluation
therefore treats refinement as a proposal:

1. Flag parents with at least 90% estimated completeness and at least 5%
   contamination.
2. Estimate a centre value
   `k0 = ceil(1 + contamination / completeness)`.
3. Evaluate `k0 - 1`, `k0`, and `k0 + 1`, bounded to `k >= 2` and routinely to
   `k <= 8`.
4. Repeat each candidate with five deterministic seeds and measure pairwise
   adjusted Rand index (ARI), both across all contigs and across contigs at
   least 5 kb long.
5. Mark a proposal stable only when median ARI is at least 0.9 in both views,
   and select the seed with the highest mean agreement to the other seeds as
   the medoid realization. Candidates implying more than eight genomes are
   routed to manual `MEGABIN_REVIEW`.

A stable result is recorded as `PENDING_QC`; it does **not** automatically
replace the parent bin. Acceptance should additionally show that marker-gene
conflicts are reduced and that splits have biologically coherent independent
support, for example taxonomy or GUNC. Composition and abundance are useful
diagnostics but are not independent evidence because VAMB already used them.

Stable medoid proposals can be materialized in a separate validation directory.
The marker-validation stage first derives a CheckM lineage marker set for the
original parent, then applies that exact marker set to the parent and every
candidate child. It retains marker IDs, gene IDs, and contig locations. A
duplicated parent marker is provisionally resolved only when all of its copies
are conserved across the children and no child retains more than one copy.
This avoids comparing contamination estimates based on different child marker
sets.

Even complete marker segregation advances a proposal only to
`PENDING_COMPLEMENTARY_QC`. Explicit child-quality thresholds, GUNC/taxonomic
coherence, controls, and a final reviewed merge remain development work. The
older forced-refinement scripts remain here for reproducibility and
experimentation, but they are not the recommended automatic production path.

## Script map

| File | Purpose | Resource |
| --- | --- | --- |
| `hierarchical_vamb_clustering_rocm.py` | Partition, local clustering, representative collection, global reconciliation, and expansion | GPU for `partition`, `local`, and `cluster-global`; CPU otherwise |
| `hierarchical_vamb.slurm` | Partition or run one local-array task | GPU |
| `hierarchical_vamb_collect.slurm` | Collect local representatives | CPU |
| `hierarchical_vamb_reconcile_cpu.slurm` | Prepare candidates, expand assignments, and export FASTAs | CPU |
| `hierarchical_vamb_reconcile_gpu.slurm` | Reconcile multi-contig representatives | GPU |
| `hierarchical_vamb_submit_sample.slurm` | Submit one complete sample as a dependency graph | CPU submitter |
| `hierarchical_vamb_submit_reconcile_checkm.slurm` | Submit reconciliation, export, and CheckM after collection | CPU submitter |
| `hierarchical_vamb_submit_wave.slurm` | Project-specific example for releasing several samples in waves | CPU submitter |
| `hierarchical_vamb_checkm.slurm` | Run CheckM with node-local intermediates | CPU/high-memory |
| `hierarchical_vamb_checkm_batch_prepare.slurm` | Create manifests for batched CheckM | CPU |
| `hierarchical_vamb_checkm_batch.slurm` | Run one CheckM batch per array task | CPU/high-memory |
| `hierarchical_vamb_checkm_batch_collect.slurm` | Validate and merge batch QA tables | CPU |
| `hierarchical_vamb_checkm_subbatch_prepare.slurm` | Extract one failed batch for smaller targeted retries | CPU |
| `hierarchical_vamb_checkm_subbatch_collect.slurm` | Reconstruct a parent QA batch from successful subbatches | CPU |
| `evaluate_hierarchical_vamb_refinement.py` | Generate and evaluate multi-seed, multi-k split proposals | CPU/GPU by subcommand |
| `hierarchical_vamb_gated_refine_cpu.slurm` | Prepare or collect gated proposals | CPU |
| `hierarchical_vamb_gated_refine_gpu.slurm` | Evaluate one candidate parent per GPU-array task | GPU |
| `hierarchical_vamb_submit_gated_refinement.slurm` | Submit gated evaluation after CheckM | CPU submitter |
| `hierarchical_vamb_submit_gated_refinement_after_prepare.slurm` | Size and submit the gated GPU array after candidate preparation | CPU submitter |
| `validate_hierarchical_vamb_proposals.py` | Materialize stable medoid proposals and report parent-matched marker conflicts | CPU |
| `hierarchical_vamb_marker_materialize.slurm` | Write stable candidate children without altering canonical bins | CPU |
| `hierarchical_vamb_marker_validate.slurm` | Run parent-matched CheckM and retain marker IDs/copy locations | CPU/high-memory |
| `refine_hierarchical_vamb.py` and `hierarchical_vamb_refine_*.slurm` | Legacy forced-refinement implementation | Experimental |
| `hierarchical_vamb_submit_refinement.slurm` | Legacy forced-refinement submitter | Experimental |

## Running one sample

The current launchers expect this layout relative to the directory from which
`sbatch` is called:

```text
vamb_crass/
  contigs.fna.gz
  SAMPLE/
    vamb/
      composition.npz
      abundance.npz
      latent.npz
```

Copy or symlink the Python and Slurm files from this directory into that
analysis directory. Then submit the dependency graph with any positive shard
count:

```bash
sbatch hierarchical_vamb_submit_sample.slurm SAMPLE 256
```

The lightweight submitter returns immediately after creating this dependency
chain:

```text
partition -> local array -> collect -> prepare -> global -> expand -> FASTA -> CheckM
```

The default workflow stops after CheckM. To also create stability-gated
refinement proposals after CheckM:

```bash
sbatch --export=ALL,HIERARCHICAL_SUBMIT_GATED_EVALUATION=1 \
  hierarchical_vamb_submit_sample.slurm SAMPLE 256
```

To submit stages manually, capture job IDs with `sbatch --parsable` and connect
deterministic stages with `--dependency=afterok:JOBID`. Do not run the local
tasks directly on a login node.

## Adapting the launchers to another system

The algorithms are dataset-independent, but the supplied Slurm wrappers record
the Pawsey environment in which they were developed. Before using them
elsewhere, review:

- Slurm account, partition, time, CPU, memory, and GPU directives;
- module names and the VAMB/CheckM environment paths;
- the `vamb_crass/SAMPLE/vamb` input convention;
- the contig FASTA path and minimum exported bin size;
- local scratch variables and filesystem quota policy.

Use scheduler allocations to determine concurrency. The workflow needs one GPU
per active partition/local/global task, not a particular total number of GPUs.
CPU-only collection, export, CheckM, and submitter jobs should remain on CPU
partitions.

## Current status

The implemented production path covers checkpoint validation, deterministic
sharding, restartable local clustering, representative collection, global
reconciliation, complete contig assignment expansion, FASTA export, and both
single-job and batched CheckM. It has been exercised on datasets large enough
to expose walltime and filesystem-file-count failure modes.

The statistical refinement layer now generates conservative, reproducible
split proposals, rejects solutions that are unstable in either contig-length
view, selects a medoid seed, and tests duplicated-marker segregation under the
same parent-derived marker set. It remains a QA and review layer: complementary
biological validation and explicit acceptance logic are not yet complete.
