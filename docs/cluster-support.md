# Adding support for another cluster

Atavide lite keeps separate scheduler scripts for each supported cluster because scheduler policy, storage, software environments, and resource limits differ in ways that are difficult to hide behind one generic script. To request a port, open a **New cluster support** issue and provide enough information to adapt one of these starting points:

- `deepthought_shortread` or `pawsey_shortread` for paired-end reads.
- `deepthought_minion` or `pawsey_minion` for single-end or long reads.

The request needs more than the scheduler name. It must identify suitable partitions or queues for ordinary CPU work, memory-heavy work, long assemblies, GPU VAMB jobs, and data transfers. It must also document resource limits, account and QoS rules, array/dependency syntax, shared and node-local temporary storage, environment setup, and a tested batch header. The issue form prompts for each of these details.

## Worked example: Pawsey Setonix

Pawsey's Setonix illustrates why a dedicated profile is useful. Setonix uses Slurm and shares compute nodes by default. Pawsey recommends explicitly requesting nodes, tasks, CPUs per task, and wall time. Shared-node jobs should request total memory per node with `--mem`; jobs requiring a complete node should use `--exclusive`.

The general-purpose partitions relevant to an atavide lite port are:

| Partition | Intended use | Maximum wall time | Resources per node | Important limits |
| --- | --- | --- | --- | --- |
| `work` | CPU production | 24 hours | 128 CPU cores, 230 GiB usable memory | General production queue |
| `long` | Long CPU production | 96 hours | 128 CPU cores, 230 GiB usable memory | One node per job; four concurrent jobs and 96 submitted jobs per user |
| `highmem` | Memory-heavy CPU production | 96 hours | 128 CPU cores, 980 GiB usable memory | One node per job; two concurrent jobs and 96 submitted jobs per user |
| `gpu` | GPU production | 24 hours | 64 CPU cores, 230 GiB usable memory, eight logical GPUs | Uses the project account with a `-gpu` suffix |
| `gpu-highmem` | GPU work needing more host memory | 48 hours | 64 CPU cores, 460 GiB usable memory, eight logical GPUs | Uses the project account with a `-gpu` suffix |
| `copy` | Large data transfers | 48 hours | 32 CPU cores, 115 GiB usable memory | Four concurrent jobs and 500 submitted jobs per user |
| `debug` | CPU development and debugging | 1 hour | 128 CPU cores, 230 GiB usable memory | Four nodes per job; one concurrent and four submitted jobs per user |
| `gpu-dev` | GPU development and debugging | 4 hours | 64 CPU cores, 230 GiB usable memory, eight logical GPUs | Two nodes per job; one concurrent and four submitted jobs per user |

These values are an example snapshot, not hard-coded pipeline assumptions. Check Pawsey's [Running Jobs on Setonix](https://pawsey.atlassian.net/wiki/spaces/US/pages/51929058/Running+Jobs+on+Setonix) and [Job Scheduling](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925964/Job+Scheduling) documentation before changing scripts, because cluster policy can change. Current limits for a user's associations can also be inspected on Setonix with:

```bash
sacctmgr show associations user=$USER cluster=setonix
```

Setonix assigns the `normal` QoS to ordinary production jobs by default. A project may use `--qos=high` for a priority boost on up to 10% of its allocation; other QoS levels are controlled by allocation state or reserved for special cases. GPU submissions use an account formed by adding `-gpu` to the base project code.

A representative shared-node CPU header is:

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

Pawsey currently recommends `--mem` rather than `--mem-per-cpu`, integer memory values, and an explicit login shell so system modules initialise correctly. Resource requests in the actual atavide lite scripts vary by pipeline stage and must remain within the chosen partition's limits.

Storage is another part of the port. Setonix profiles use shared scratch under `/scratch/$PAWSEY_PROJECT/$USER`; scratch files are subject to a 21-day purge policy. Project data that must persist is staged to or from Pawsey's Acacia object storage, with large transfers assigned to the `copy` partition. Conda environments and databases used by the Pawsey profiles are also placed under the project scratch hierarchy, so the scripts include checks and recreation steps for purged environments.

This example supplies the information needed to decide how to translate scheduler headers and paths. A complete new-cluster request should additionally supply a tested batch script, exact job-array and dependency syntax, software/module versions, GPU runtime requirements, filesystem quotas, and cluster-specific operational restrictions.
