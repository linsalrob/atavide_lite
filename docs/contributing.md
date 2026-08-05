# Contributing

Contributions are welcome from people of every background and experience level, including first-time open-source contributors. Bug reports, questions, documentation, tests, utilities, and cluster profiles all improve the project.

Communicate respectfully and assume good intent. It is fine not to know the repository, Git, an HPC scheduler, or a bioinformatics tool yet. Ask for help when something is unclear, avoid sharing sensitive data or credentials, and follow institutional policy.

## Start with an issue

Search existing issues before starting. Use the bug or feature template for general work. For a new HPC system, use the **New cluster support** form and read [Adding support for another cluster](cluster-support.md).

An issue gives maintainers a chance to confirm scope, point to related work, and identify scientific or operational constraints before implementation.

Include evidence that lets another person understand the request:

- repository commit/profile and relevant command;
- expected and observed behaviour;
- scheduler/cluster and redacted resource header;
- minimal logs or traceback;
- reproducible sample naming/layout without sequence data; and
- official documentation supporting cluster-specific values.

Never post passwords, tokens, private project identifiers, controlled data, or sensitive paths.

## Create a focused branch

Fork the repository if necessary, clone it, and create a branch from the current default branch:

```bash
git clone https://github.com/<your-user>/atavide_lite.git
cd atavide_lite
git switch -c docs/clear-topic-name
```

Use a cluster-specific branch such as `cluster/examplecluster_shortread` for a profile port. Keep unrelated changes on separate branches.

## Make changes that remain inspectable

- Preserve the visible stage boundaries and control-file contract unless the proposal intentionally changes them.
- Do not put usernames, account codes, or private filesystem paths into shared scripts.
- Keep application thread counts consistent with scheduler allocations.
- Update documentation when commands, outputs, dependencies, or profile behaviour change.
- Prefer a small, reviewable change over a large unrelated rewrite.
- Do not commit input data, generated outputs, databases, scheduler logs, or environments.

## Using agentic AI

Coding assistants can help inventory scripts and perform repetitive profile adaptation. Give the assistant the issue, official HPC documentation, and the closest existing profile. One useful request is:

> Create `examplecluster_shortread` from `pawsey_shortread`. Use the cluster-support issue and official scheduler documentation to adapt every directive, resource request, storage path, environment step, and README command. Do not modify existing profiles. List anything that cannot be verified and run local checks without submitting a real job.

Ask the assistant to inspect the complete source profile and show its diff. Then personally verify every queue, QoS, account, time, CPU, memory, GPU, path, module, and command. AI can invent plausible but invalid HPC details. Do not authorise it to submit jobs, access sensitive data, or expose credentials unless those actions are explicitly intended and understood.

## Validate locally and on the cluster

Before committing:

1. Run `git diff --check` and review the entire diff.
2. Run the repository unit tests.
3. Parse or build documentation and check links.
4. Run `bash -n` on changed shell and scheduler scripts.
5. Search new profiles for stale cluster names, paths, accounts, partitions, and environment variables.
6. Check resource requests against current official policy.
7. Submit a harmless minimal job and then a small non-sensitive pipeline example when cluster access is available.

Document what was tested and what remains unverified.

## Commit and open a pull request

Stage only intended paths:

```bash
git status --short
git add docs/ examplecluster_shortread/
git commit -m "Add ExampleCluster short-read profile"
git push -u origin cluster/examplecluster_shortread
```

Open a draft pull request through GitHub or the CLI:

```bash
gh pr create --draft --fill
```

The pull request should link its issue and explain:

- the problem and approach;
- the source profile, if any;
- scheduler, storage, environment, and scientific changes;
- documentation or evidence supporting decisions;
- checks and sample jobs completed; and
- remaining risks or unverified assumptions.

Keep the PR in draft while implementation or target-cluster testing remains. Respond to review comments with follow-up commits so the review history remains useful. Maintainers decide when the change is ready to merge.

## Documentation contributions

Install the documentation dependencies, then build with warnings treated as errors:

```bash
python -m pip install -r docs/requirements.txt
sphinx-build -W --keep-going -b html docs docs/_build/html
```

Read the rendered pages, not only the Markdown. Navigation, code blocks, tables, cross-references, and external links should remain clear on narrow screens and for readers unfamiliar with the project.
