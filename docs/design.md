# Design and rationale

## Explicit stages instead of one opaque run

The original, more comprehensive atavide pipeline used Snakemake as its main workflow manager. That approach automated scheduling effectively, but operational experience with large and heterogeneous datasets showed a cost: failures at different stages were difficult to isolate, and adapting global resource settings to different clusters was cumbersome.

`atavide_lite` therefore treats the HPC scheduler as the top-level orchestrator. Users submit visible scripts for each major stage and connect them with success dependencies such as Slurm's `afterok`. Snakemake is still used where its rule graph is helpful, including taxonomy summarisation, but it no longer hides the entire end-to-end execution.

This is an intentional trade-off. Users issue more explicit commands, but gain direct control over restart points, resource requests, logs, and cluster-specific behaviour.

## Failure containment and restartability

Every pipeline stage produces durable outputs before downstream work starts. A failed MMseqs2 array task, for example, can be investigated and resubmitted without rerunning quality control or host removal. Scheduler output and error files preserve the job and array identifiers, making failures traceable to a stage and sample.

Dependencies express ordering without requiring a long-lived process on a login node. Downstream jobs start only after the required upstream job succeeds. Where a downstream diagnostic should run even after failure, individual profiles may deliberately use `afterany`; users should inspect those choices before submission.

## Separate profiles instead of universal scheduler abstraction

The project maintains distinct directories because cluster differences extend beyond translating `#SBATCH` to `#PBS`. Profiles encode:

- paired versus single-read filename handling;
- Slurm or PBS submission and dependency syntax;
- partition, queue, QoS, account, and GPU directives;
- job-array limits;
- shared, node-local, and persistent storage paths;
- scratch staging and cleanup;
- shell, modules, Conda/Mamba, and GPU runtime setup; and
- resource estimates for particular pipeline stages.

A universal template would either expose all of these differences as configuration or conceal them behind generated code. Separate readable profiles make local policy visible and allow maintainers to optimise the systems they actually operate.

## A small analysis contract

The profiles share a lightweight contract rather than a large configuration framework:

- `DEFINITIONS.sh` records dataset-specific paths and naming decisions.
- `R1_reads.txt` controls paired-end samples; `reads.txt` controls single/long-read samples.
- Scheduler array indices select entries from those files.
- Stable directory names connect outputs from one stage to inputs of the next.
- Job IDs returned by the scheduler create dependency chains.

These plain-text control files become part of the analysis record and are easy to inspect, version, or regenerate.

## Similar processing across sequencing technologies

Paired short reads and long reads require different filename handling and may use different quality thresholds, but their downstream questions are similar. The project therefore keeps stage names and overall logic aligned where practical: quality control, optional host removal, conversion or ORF preparation, taxonomy, function, assembly, and binning.

The separation is retained where it improves clarity. Paired reads are tracked as R1/R2 mates, while single or long reads use one sample file. Defaults should be reviewed for the technology and biological question; a convenient profile name is not a substitute for checking tool parameters.

## Read-based and assembly-based branches

Read classification and assembly answer different questions, so neither is treated as a replacement for the other. Read-level MMseqs2 results support abundance summaries. Assembly and VAMB support contig- and genome-level exploration. Because the branches diverge after cleaned reads are available, assembly work does not need to block read-based annotation unless a specific downstream dependency requires it.

## Resource-aware data movement

Large databases and temporary files make storage architecture part of the computation. Deepthought profiles can use fast node-local `$BGFS` space, while Pawsey profiles use shared `/scratch` and account for its purge policy. NCI jobs can use `$PBS_JOBFS`. Moving temporary I/O to suitable storage can materially improve runtime, but every profile must also copy required results back to durable, shared storage.

## Conservative use of automation and AI

Repetition makes profile adaptation suitable for coding assistants, but scheduler directives can be plausible and still wrong. The contribution process requires official cluster documentation, a tested minimal batch header, diff review, and human verification of every queue, account, resource, path, and module choice. Automation assists the port; it does not validate institutional policy.
