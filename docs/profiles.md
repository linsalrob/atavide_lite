# Choosing a profile

Select a profile using two independent questions: how are the reads represented, and which cluster's scheduler/storage model is closest?

## Read type

| Input | Control file | Closest profiles |
| --- | --- | --- |
| Paired-end short reads (`R1` and `R2`) | `R1_reads.txt` | `pawsey_shortread`, `deepthought_shortread`, `nci_pbs/fasta` |
| Single-end or long reads | `reads.txt` | `pawsey_minion`, `deepthought_minion` |

Illumina single-end data uses the single-read control model even though the reads are short. Conversely, paired data downloaded with `_1.fastq.gz` and `_2.fastq.gz` may need renaming to `_R1.fastq.gz` and `_R2.fastq.gz` before using a paired profile.

## Supplied cluster profiles

### Pawsey Setonix

`pawsey_shortread/` and `pawsey_minion/` use Slurm and shared project scratch under `/scratch/$PAWSEY_PROJECT/$USER`. They include environment checks because files on scratch are subject to purge. Large persistent data should be staged to and from Pawsey's Acacia object storage, and transfer work should follow Pawsey's current partition policy.

Choose these profiles when using Setonix. They may also be the closest starting point for a different Slurm cluster with shared scratch, but they are not generic: project variables, account rules, partitions, GPU modules, paths, and purge assumptions must be reviewed.

### Flinders Deepthought

`deepthought_shortread/` and `deepthought_minion/` use Slurm. Some stages take advantage of fast node-local storage exposed through `$BGFS`. That storage is not a shared hand-off location, so scripts stage inputs in and copy required results back out.

Choose these profiles on Deepthought or as the starting point for a cluster whose node-local temporary storage is central to performance.

### NCI Gadi

`nci_pbs/fasta/` contains PBS scripts for Gadi. It reflects PBS directives, project charging, queue selection, Gadi array constraints, shared `/g/data`, and compute-node-local `$PBS_JOBFS`.

The NCI profile is less complete than the paired Slurm profiles and should be treated as an existing PBS implementation to inspect, not as a scheduler-neutral conversion layer.

### Generic and component workflows

`bash/atavide_lite.sh` demonstrates the main commands in one Bash script. It is useful for learning and adaptation, but the repository explicitly recommends running and checking stages one at a time.

`summarise_taxonomy/` and `minion_taxonomy_and_function/` are component workflows. Use them when only those functions are needed or when integrating atavide lite outputs into a larger analysis.

## Selection checklist

Before running a profile, confirm all of the following:

- filenames match the paired or single-read assumptions;
- scheduler directives are valid on the target cluster;
- account, project, partition/queue, and QoS values are permitted;
- requested time, CPUs, memory, nodes, GPUs, and array sizes are within limits;
- all hard-coded and environment-derived paths exist on compute nodes;
- temporary output is copied back before node-local storage disappears;
- the job shell initialises Conda and modules correctly;
- database versions and locations are correct; and
- profile README commands reference scripts that exist in that directory.

If no supplied profile passes this checklist, follow [Adding support for another cluster](cluster-support.md) instead of modifying a maintained profile in place.
