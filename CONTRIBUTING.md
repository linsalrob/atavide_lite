# Contributing to atavide lite

Thank you for helping improve atavide lite. Contributions are welcome from people of every background and experience level, including first-time open-source contributors. Clear questions, documentation fixes, bug reports, tests, and cluster profiles are all valuable contributions.

Please communicate respectfully and assume good intent. It is fine not to know the repository, Git, a scheduler, or a bioinformatics tool yet. Ask questions when something is unclear, avoid sharing sensitive research data or credentials, and follow your institution's policies when describing its systems.

## Start with an issue

Before changing code, search the existing issues for related work. For bugs and general improvements, use the relevant issue template. For support for a new HPC system, open the **New cluster support** form and give as much detail as possible about:

- Slurm, PBS, or the scheduler your cluster uses.
- CPU, high-memory, GPU, long-running, debug, and transfer partitions or queues.
- Maximum wall times, CPUs, memory, nodes, GPUs, job arrays, and queue limits.
- Required project/account, QoS, queue, GPU, export, or other submission directives.
- Shared and node-local scratch storage, long-term storage, quotas, and purge policies.
- Modules, Conda/Mamba, containers, GPU runtimes, and other environment setup.
- Submission, dependency, array, status, accounting, and cancellation commands.

Link to current cluster documentation and include a small batch script that has run successfully. Ask your HPC support team for missing details. Do not include passwords, tokens, private project identifiers, sensitive paths, or research data. The [cluster support guide](docs/cluster-support.md) uses Pawsey Setonix as an example of a complete system description.

Starting with an issue lets maintainers and other users confirm the intended profile, identify policy constraints, and avoid duplicating work.

## Make a branch

Fork the repository if you do not have write access, clone your fork, and create a focused branch from the current default branch:

```bash
git clone https://github.com/<your-user>/atavide_lite.git
cd atavide_lite
git switch -c cluster/<short-cluster-name>
```

Keep unrelated changes on separate branches. If the work takes time, regularly incorporate upstream changes using the Git workflow you are comfortable with.

## Adapt the closest profile

Choose a starting directory by input type and cluster behaviour:

- Use `deepthought_shortread` or `pawsey_shortread` for paired-end reads.
- Use `deepthought_minion` or `pawsey_minion` for single-end or long reads.
- Prefer the Pawsey profiles for systems with shared scratch similar to `/scratch`.
- Prefer the Deepthought profiles when node-local temporary storage is the closer model.

Copy the entire closest directory to a clearly named new directory, then update the copy. Do not change an existing cluster profile merely to make a different cluster work. For example:

```bash
cp -R pawsey_shortread examplecluster_shortread
```

Review every scheduler script, definition file, environment check, path, and README in the new directory. Translate scheduler directives, resource requests, arrays, dependencies, output paths, storage locations, software activation, modules, GPU setup, and submission examples. Keep the pipeline stages consistent with the source profile unless the issue describes a deliberate functional change.

## Using agentic AI

An agentic coding assistant can perform repetitive copying and scheduler-header changes, but it needs precise source material and human review. Give the assistant the cluster-support issue, links to official HPC documentation, and the name of the closest profile. A useful request is:

> Create `examplecluster_shortread` from `pawsey_shortread`. Use the attached cluster-support issue and official scheduler documentation to adapt every scheduler directive, resource request, storage path, environment setup step, and README command. Do not modify existing cluster profiles. List any policy or parameter you cannot verify, and run repository checks without submitting a real cluster job.

Ask the assistant to inspect the whole source directory before editing and to show a diff afterward. Then personally verify every generated resource request and command against current cluster documentation. AI output can contain plausible but invalid queue names, directives, paths, module versions, or limits. Never let an assistant submit jobs, expose credentials, or access sensitive data unless that action is explicitly authorised and understood. Test with a harmless minimal job before using real data or expensive allocations.

## Check the contribution

Before committing:

1. Search the new directory for old cluster names, account codes, paths, partitions, and environment variables.
2. Check that each requested wall time, CPU count, memory value, node count, GPU request, and array size is allowed in its selected queue.
3. Check shell syntax and run the repository's relevant automated tests.
4. Confirm that README commands use the new directory and scheduler.
5. Run one minimal job and then a small, non-sensitive pipeline example on the target cluster when access is available.
6. Review `git diff` and make sure the branch contains no generated data, scheduler output, credentials, or unrelated files.

If cluster access is unavailable, say exactly what was checked locally and what still needs testing. That does not prevent early discussion or a draft pull request.

## Commit and open a pull request

Stage only the intended files, write a concise commit message, and push your branch:

```bash
git status --short
git add examplecluster_shortread
git commit -m "Add ExampleCluster short-read profile"
git push -u origin cluster/examplecluster
```

Open a pull request against the upstream repository using GitHub's website or GitHub CLI:

```bash
gh pr create --draft --fill
```

In the pull request, link the cluster-support issue and explain:

- Which existing profile was used as the starting point.
- The scheduler, queues, storage, and software-environment changes made.
- Which documentation and cluster policies support the chosen values.
- Tests and sample jobs that passed.
- Anything that remains unverified or needs a maintainer or HPC administrator to review.

Use a draft pull request while work or cluster testing remains. Respond to review comments with follow-up commits rather than replacing discussion history. A maintainer can help decide when the contribution is ready to merge.
