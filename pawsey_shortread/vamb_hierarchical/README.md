# Hierarchical VAMB clustering

This workflow is a scalable, stand-alone continuation from a VAMB feature
checkpoint. It requires only these files from a VAMB directory:

- `composition.npz`
- `abundance.npz`
- `latent.npz`

It does not use partial cluster tables or modify the standard VAMB resume
workflow.

## Why use a hierarchy?

VAMB's `ClusterGenerator` is a destructive, greedy medoid algorithm: each
accepted bin changes the set of contigs available to the next decision. Its
global clustering loop is therefore serial. For a 22-million-contig latent
matrix, a normal ROCm-resumed run made too little progress to fit a 24-hour
walltime. Splitting matrix-vector distance calculations across GPUs did not
help: GPU-to-GPU transfers outweighed the small parallel kernels.

The hierarchy makes the expensive clustering local while preserving a
statistical reconciliation path:

1. GPU spherical k-means partitions the 32-dimensional VAMB latent space into
   256 deterministic, spatially coherent shards.
2. An unthrottled Slurm array runs normal VAMB clustering inside each shard.
3. Each local bin is represented by its latent centroid, total bases, and
   number of contigs.
4. Multi-contig local bins are globally reconciled with VAMB. Local singleton
   bins are retained explicitly, matching VAMB loner semantics rather than
   silently dropping contigs.
5. The global mapping is expanded back to a standard
   `clustername<TAB>contigname` table, then the existing
   `vamb_create_fasta.py` helper writes FASTA bins.

This is an approximation to a monolithic VAMB clustering. It is designed to
make the trade-off explicit, preserve every contig, and retain conventional
VAMB-compatible output tables.

## Scripts

| File | Purpose | Queue |
| --- | --- | --- |
| `hierarchical_vamb_clustering_rocm.py` | Partitioning, local clustering, representative collection, global reconciliation, and expansion | GPU for `partition`, `local`, and `cluster-global`; CPU otherwise |
| `hierarchical_vamb.slurm` | Feature-checkpoint partition and local-array launcher | GPU |
| `hierarchical_vamb_collect.slurm` | Collect local representatives and member counts | CPU `work` |
| `hierarchical_vamb_reconcile_cpu.slurm` | Prepare global candidates, expand memberships, and write FASTAs | CPU `work` |
| `hierarchical_vamb_reconcile_gpu.slurm` | Cluster multi-contig representatives globally | GPU |
| `hierarchical_vamb_submit_wave.slurm` | Releases sample waves after prior collection jobs complete | CPU `work` |
| `hierarchical_vamb_checkm.slurm` | CheckM with deep-tree recursion handling and tabular QA output | CPU `highmem` |
| `hierarchical_vamb_checkm_batch_*.slurm` | Manifest preparation, independent CheckM batches, and QA-table merge | CPU `work` / `highmem` |
| `refine_hierarchical_vamb.py` | CheckM-guided targeted refinement of contaminated bins | CPU/GPU by stage |
| `hierarchical_vamb_refine_cpu.slurm` | Extract candidates, merge refined assignments, and write refined FASTAs | CPU `work` |
| `hierarchical_vamb_refine_gpu.slurm` | Re-cluster one contaminated parent bin per GPU-array task | GPU |
| `hierarchical_vamb_submit_refinement.slurm` | Derives the refinement array size from CheckM QA and submits one pass | CPU `work` |

## CheckM-guided refinement

The global representative reconciliation is deliberately permissive enough to
restore bins split at shard boundaries. In some cases, that can over-merge
several genomes. The post-CheckM refinement pass corrects only this failure
mode without changing unaffected assignments:

1. CheckM writes `checkm/qa.tsv` in tabular form.
2. Bins with completeness at least 90% and contamination at least 5% are
   selected as likely merged bins.
3. Their original contigs are extracted from the complete assignment table and
   force-partitioned by GPU spherical k-means in the original VAMB embedding.
   CheckM's marker-copy burden calibrates the number of groups: approximately
   `1 + contamination / completeness`, bounded between 2 and 64.
4. The resulting sub-bins replace only those parent assignments in
   `refinement/refined_unsplit.tsv`; every other contig is retained unchanged.
5. Refined FASTAs are checked with CheckM again.

The automatic submitter runs one GPU array task per selected parent bin and
does not impose an artificial array throttle. If no bins meet the criteria, it
completes without launching a GPU array.

CheckM runs use node-local temporary storage for their per-bin HMM
intermediates and retain only `qa.tsv` in the project directory. This avoids
exhausting the Lustre per-user file-count quota for samples with thousands of
bins.

For exceptionally large bin sets such as Whaleshark, use the batch launchers
to validate independent groups of roughly 1,000 bins per high-memory array
task. Their QA tables are merged into the same `checkm/qa.tsv` consumed by the
refinement workflow, so batching changes the scheduling and failure domain,
not the downstream interface.

## Example commands

Use a completed VAMB directory at `vamb_crass/SAMPLE/vamb`.

```bash
# 1. Partition the latent space.
sbatch --partition=gpu-dev --gres=gpu:1 --time=04:00:00 \
  --export=ALL,HIERARCHICAL_STAGE=partition,HIERARCHICAL_SHARDS=256 \
  hierarchical_vamb.slurm SAMPLE

# 2. After partitioning, cluster every shard. Do not add an artificial %N
# array throttle; Slurm/account policy controls concurrency.
sbatch --array=0-255 --partition=gpu --gres=gpu:1 --time=1-00:00:00 \
  --export=ALL,HIERARCHICAL_STAGE=local,HIERARCHICAL_SHARDS=256 \
  hierarchical_vamb.slurm SAMPLE

# 3. Collect representatives without consuming a GPU.
sbatch hierarchical_vamb_collect.slurm SAMPLE

# 4. Reconcile multi-contig representatives, expand the final table, and
# write FASTA bins.
sbatch hierarchical_vamb_reconcile_cpu.slurm SAMPLE prepare
sbatch hierarchical_vamb_reconcile_gpu.slurm SAMPLE
sbatch hierarchical_vamb_reconcile_cpu.slurm SAMPLE expand
sbatch hierarchical_vamb_reconcile_cpu.slurm SAMPLE fasta
```

For production, connect these commands with `--dependency=afterok:JOBID`.

## Whaleshark example results

The initial Whaleshark run used 22,193,255 contigs and completed the first
two hierarchy stages successfully:

| Measurement | Result |
| --- | ---: |
| Latent-space shards | 256 |
| Shard size (min / median / max) | 32,051 / 84,817 / 168,297 contigs |
| Partition walltime | 48m 57s on one MI250X GPU |
| Local array result | 256 / 256 shards completed; no failures |
| Local representatives | 15,351,451 |
| Multi-contig representatives | 238,852 (1.56%) |
| Singleton representatives | 15,112,599 |
| Representative collection walltime | 12m 52s on CPU |

The global multi-contig reconciliation and final FASTA-generation stages are
submitted as dependent jobs. Final FASTA output follows the established VAMB
pattern:

```bash
python ~/atavide_lite/bin/vamb_create_fasta.py \
  vamb_crass/contigs.fna.gz \
  vamb_crass/Whaleshark/hierarchical_vamb/hierarchical_unsplit.tsv \
  20000 \
  vamb_crass/Whaleshark/hierarchical_vamb/hierarchical_unsplit_bins
```
