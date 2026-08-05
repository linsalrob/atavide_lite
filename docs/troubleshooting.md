# Troubleshooting

## A practical diagnosis order

When a job fails:

1. Identify the job and array task in the scheduler.
2. Read the complete `.err` and `.out` files.
3. Find the last command printed or the line reported by the shell error trap.
4. Check whether expected inputs exist and are non-empty.
5. Confirm `DEFINITIONS.sh`, the sample list entry, and derived filename.
6. Check scheduler exit state, elapsed time, maximum memory, and node failure information.
7. Reproduce with one small sample or a harmless diagnostic job.
8. Fix the cause and resubmit only the failed stage and its blocked dependants.

Do not immediately add memory or time. Incorrect paths, filenames, databases, or environment activation are more common and will fail again with a larger allocation.

## Unbound variables

Scripts use `set -u`, so a missing variable stops the job. Verify that:

- `DEFINITIONS.sh` exists in the submission working directory;
- assignments have no spaces around `=`;
- required values such as `FILEEND` and `HOSTREMOVED` are present;
- cluster variables such as `PAWSEY_PROJECT`, `BGFS`, or `PBS_JOBFS` exist in batch jobs; and
- the job was submitted from the intended analysis directory.

## No samples or wrong array indices

Print the control file with line numbers:

```bash
nl -ba R1_reads.txt
```

Slurm arrays in the supplied scripts are normally one-based because selection uses `head -n $SLURM_ARRAY_TASK_ID`. An array starting at zero can select no record. Remove blank lines and regenerate the list with `sort` for deterministic ordering.

## Missing paired reads

Paired scripts often replace `_R1` with `_R2`. Check capitalization, lane/run suffixes, and whether `FILEEND` removes the correct suffix. Confirm each R1-derived R2 path exists before submitting an array.

## Conda activation fails in a batch job

An environment working interactively does not prove the batch shell initialises it. Inspect the script's shebang and Conda hook. Confirm the environment path is visible from compute nodes and has not been purged. On Pawsey, run the supplied environment check/recreation workflow.

Avoid editing global shell startup files as a first response; make environment setup explicit and reproducible inside the profile.

## Tool not found or wrong version

Print versions after activation in a small scheduler job:

```bash
fastp --version
minimap2 --version
samtools --version
mmseqs version
```

Check for a site module overriding the Conda executable and verify GPU/runtime modules match the VAMB build.

## Job exceeded wall time

Determine whether the job made progress and whether its runtime scales with reads, database size, or assembly complexity. Use scheduler accounting from successful small jobs to estimate a new request. If the required time exceeds the current partition limit, choose an authorised long partition and adapt its node/memory constraints rather than requesting an invalid duration.

## Out of memory

Confirm the application actually received and used the requested CPUs and memory. MMseqs2, MEGAHIT, mapping, and cross-sample VAMB preparation can scale strongly with input size. Possible responses include reducing concurrency, splitting work, using per-sample assembly, selecting an authorised high-memory queue, moving temporary data to suitable storage, or increasing memory based on measured peak use.

Do not assume `--mem-per-cpu` and `--mem` are interchangeable; follow cluster policy.

## Files disappeared from scratch

Scratch is temporary by design. Check the provider's purge policy and object/persistent storage. Recreate environments or databases from recorded definitions, then stage required inputs back. Build transfer and archival jobs into the workflow rather than relying on manual access to reset file age.

## MMseqs2 database not found

Check the expanded database path from inside a compute job, not only on the login node. Confirm all files that make up the MMseqs2 database are present and readable, and that the script's UniRef name matches the installed database.

## TaxonKit errors or incomplete taxonomy

Verify `TAXONKIT_DB` points to a directory containing `names.dmp`, `nodes.dmp`, `delnodes.dmp`, and `merged.dmp`. Record the dump date. Use a PyTaxonKit version compatible with current NCBI taxonomy formatting.

## Empty output after fastp

Read the fastp JSON/HTML report. The minimum-length or maximum-N threshold may exclude an older or unusual library. Change thresholds only after deciding they remain scientifically appropriate for downstream tools, and record the change.

## Dependency will never run

Inspect the prerequisite job's terminal state. A job submitted with `afterok` remains blocked or is cancelled when the prerequisite fails. After fixing and resubmitting the failed job, submit a new downstream dependency using the new job ID; changing files does not alter an already-recorded dependency.

## Asking for help

Open a bug report with the repository commit, profile, scheduler/cluster, exact submission command, relevant script header, redacted `DEFINITIONS.sh`, sample naming pattern, scheduler state/accounting, and the smallest relevant log excerpt. Remove credentials, private project codes, sensitive paths, and sequence data.
