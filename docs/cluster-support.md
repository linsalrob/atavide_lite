# Adding support for another cluster

Atavide lite keeps separate scheduler profiles because resource policy, storage, software setup, and operational limits differ even between systems using the same scheduler. A successful port preserves the scientific stage contract while translating those system-specific details explicitly.

## Begin with a cluster-support issue

Search existing issues, then open the repository's **New cluster support** form. A complete request should contain:

- cluster and organisation name;
- links to current official scheduler, storage, and software documentation;
- scheduler family and exact version when relevant;
- sequencing input model to support;
- all CPU, high-memory, GPU, long, debug, and transfer partitions/queues;
- maximum wall time, cores, usable memory, nodes, GPUs, arrays, running jobs, and submitted jobs;
- account/project, QoS, partition/queue, constraint, reservation, export, and GPU rules;
- shared scratch, node-local scratch, persistent storage, quotas, and purge periods;
- shell, module, Conda/Mamba, container, and GPU-runtime setup;
- submit, parsable job ID, array, dependency, status, accounting, cancellation, and task-launch syntax; and
- a small batch script known to run successfully.

Ask HPC support for anything uncertain. Remove credentials, private account codes, sensitive paths, and research data.

## Choose the closest source profile

| Need | Starting points |
| --- | --- |
| Paired-end reads, shared scratch, Slurm | `pawsey_shortread` |
| Single/long reads, shared scratch, Slurm | `pawsey_minion` |
| Paired-end reads, node-local scratch, Slurm | `deepthought_shortread` |
| Single/long reads, node-local scratch, Slurm | `deepthought_minion` |
| Paired reads, PBS | `nci_pbs/fasta` |

Copy the entire closest directory to a new, clearly named profile. Do not change an existing maintained profile merely to support a different system.

```bash
cp -R pawsey_shortread examplecluster_shortread
```

## Porting checklist

Review every file in the copy. At minimum translate:

1. **Batch headers:** shebang, account/project, queue/partition, QoS, time, tasks, CPUs, memory, nodes, GPUs, arrays, output/error, exports, constraints, and exclusivity.
2. **Submission examples:** job-ID parsing, arrays, dependency conditions, throttling, and cancellation.
3. **Storage:** analysis paths, databases, shared and node-local scratch, copying, cleanup, persistent storage, and purge policy.
4. **Environment:** login shell, modules, Conda hooks/prefixes, containers, compilers, and GPU runtime.
5. **Resources:** make application thread counts agree with scheduler allocations and map memory-heavy, long, GPU, and transfer stages to valid queues.
6. **Control files:** retain the paired `R1_reads.txt` or single `reads.txt` contract and appropriate `DEFINITIONS.sh` examples.
7. **Documentation:** update every path, variable, command, queue, environment, and profile name in the copied README.

Search for stale source-cluster terms after editing:

```bash
rg -n 'PAWSEY|Pawsey|BGFS|PBS_JOBFS|/scratch|--partition|--qos|#PBS' \
  examplecluster_shortread
```

The presence of a match is not automatically wrong; each match must be understood.

## Resource mapping by stage

A cluster description is not complete until the main workload classes have destinations:

| Workload | Questions to answer |
| --- | --- |
| fastp and light summaries | Ordinary CPU queue? Array/concurrency limit? |
| Host mapping and read mapping | CPU/thread count, memory, local I/O strategy? |
| MMseqs2 | High memory or standard memory? Temporary storage? Maximum duration? |
| MEGAHIT | Per-sample versus cross-assembly, memory, restart, long queue? |
| VAMB | CPU or GPU build? GPU type/count/runtime/account? |
| Database download and archival | Transfer nodes/queue? Outbound network? Persistent target? |

Measure small representative jobs and revise defaults based on scheduler accounting. Requesting the maximum available resource for every stage wastes allocations and may delay scheduling.

## Worked example: Pawsey Setonix

Setonix demonstrates why a separate profile is needed. It uses Slurm, shares compute nodes by default, and recommends explicit node, task, CPU, and wall-time requests. Shared-node jobs request total node memory with `--mem`; complete-node jobs use `--exclusive`.

General-purpose partitions relevant to atavide lite include:

| Partition | Use | Maximum wall time | Resources per node | Selected limits |
| --- | --- | --- | --- | --- |
| `work` | CPU production | 24 h | 128 cores, 230 GiB usable memory | General production |
| `long` | Long CPU work | 96 h | 128 cores, 230 GiB | One node/job; 4 running and 96 submitted/user |
| `highmem` | Memory-heavy CPU work | 96 h | 128 cores, 980 GiB | One node/job; 2 running and 96 submitted/user |
| `gpu` | GPU production | 24 h | 64 cores, 230 GiB, 8 logical GPUs | GPU project account required |
| `gpu-highmem` | GPU with more host memory | 48 h | 64 cores, 460 GiB, 8 logical GPUs | GPU project account required |
| `copy` | Large transfers | 48 h | 32 cores, 115 GiB | 4 running and 500 submitted/user |
| `debug` | CPU development | 1 h | 128 cores, 230 GiB | 4 nodes/job; 1 running and 4 submitted/user |
| `gpu-dev` | GPU development | 4 h | 64 cores, 230 GiB, 8 logical GPUs | 2 nodes/job; 1 running and 4 submitted/user |

These limits are a documentation snapshot, not pipeline constants. Confirm them using Pawsey's [Running Jobs on Setonix](https://pawsey.atlassian.net/wiki/spaces/US/pages/51929058/Running+Jobs+on+Setonix), [Job Scheduling](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925964/Job+Scheduling), and current association data:

```bash
sacctmgr show associations user="$USER" cluster=setonix
```

Ordinary production uses the default `normal` QoS. Pawsey documents `--qos=high` as a limited priority boost. GPU jobs charge an account formed by adding `-gpu` to the base project code.

A representative shared-node header is:

```bash
#!/bin/bash --login
#SBATCH --account=<project>
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=58880M
#SBATCH --time=01:00:00
```

Pawsey recommends integer memory values and currently favours `--mem` over `--mem-per-cpu`. Shared project scratch is under `/scratch/$PAWSEY_PROJECT/$USER` and is subject to a 21-day purge policy. Persistent data are staged through Acacia object storage, with large transfers assigned to `copy`. Environments and databases stored on scratch must be reproducibly recreatable.

## Validate the new profile

Validation should progress from harmless to representative:

1. Check shell syntax for every batch script.
2. Check paths and stale identifiers with `rg`.
3. Submit a minimal job that prints hostname, storage, environment, and versions.
4. Submit one array task and one dependency chain.
5. Process a small non-sensitive sample through every stage.
6. Record scheduler time, peak memory, CPU efficiency, and GPU use.
7. Test restart behaviour after an intentionally stopped disposable job.
8. Ask a local HPC administrator or experienced user to review directives.

If cluster access is unavailable, state that clearly in the pull request and open it as a draft.
