# Pipeline stages

The pipeline contains a shared preparation phase followed by read-based and assembly-based branches. Script names vary slightly by profile; inspect the selected directory and its README before copying a command.

## Dependency model

For Slurm, `sbatch --parsable` returns a job ID that can be used in later submissions:

```bash
FIRST_JOB=$(sbatch --parsable first.slurm)
SECOND_JOB=$(sbatch --parsable \
  --dependency="afterok:${FIRST_JOB}" second.slurm)
```

`afterok` starts the dependent job only if the prerequisite exits successfully. Array semantics and PBS dependency syntax differ. A dependency cannot verify scientific validity: inspect outputs before continuing when first using a stage or dataset.

## 1. Quality control with fastp

The fastp stage removes adapters and filters reads. Current profile defaults differ by input type. The paired short-read scripts generally allow at most one uncalled base and require a minimum length near 100 bp; single/ONT examples may allow two uncalled bases and use a shorter minimum in some profiles.

These are project defaults, not universal biological rules. Short archived reads can all fail a 100 bp threshold. Review the exact `fastp` command, adapter FASTA, JSON/HTML report, retained-read counts, and length distribution before accepting the output.

Typical outputs include cleaned FASTQ files under `fastq_fastp/` and reports under `fastp_output/`.

## 2. Optional host-read separation

`host_removal` uses minimap2 to align cleaned reads to `HOSTFILE` and samtools flags to separate mapped and unmapped primary reads that pass platform/vendor checks and are not marked as PCR/optical duplicates. Paired profiles additionally retain mate identity when writing R1 and R2 outputs.

The reference can contain any sequence set that should be removed from downstream analysis. The mapped and unmapped reads are stored separately; downstream stages consume `HOSTREMOVED`.

Validate:

- reference identity, version, and indexing;
- mapping and retained-read percentages;
- mate pairing after filtering;
- whether unexpectedly high host content indicates a sample or reference problem; and
- data-governance requirements for retained host reads.

To skip host alignment, configure downstream input to use the fastp output as described by the profile rather than inventing a dummy reference.

## 3. FASTQ-to-FASTA conversion

MMseqs2 `easy-taxonomy` consumes FASTA input. The repository's compiled `fastq2fasta` utility converts cleaned reads efficiently and can process a directory. Paired profiles preserve `/1` and `/2` identity in headers where required.

Confirm that the binary was compiled, output files are non-empty, compression is as expected, and derived sample names match `FILEEND`/`FAFILEEND`.

## 4. Read-based taxonomy with MMseqs2

MMseqs2 `easy-taxonomy` compares reads with a UniRef database and assigns taxonomy using a lowest-common-ancestor approach after its prefiltering and alignment stages. This is commonly the most computationally demanding read-based stage because both query volume and database size are large.

Current Pawsey scripts provide UniRef50 and UniRef100 variants, while the current project methods favour UniRef100 for greater coverage. The database choice affects runtime, storage, hit rate, identifiers, and interpretation. Do not mix results from different UniRef releases without recording and accounting for the difference.

MMseqs2 uses substantial temporary storage. Select a long or high-memory queue only when its limits and charging model fit the job, and place temporary files on storage designed for heavy I/O.

## 5. Taxonomy completion and summaries

The `summarise_taxonomy` component uses TaxonKit through PyTaxonKit to expand taxonomic lineages from MMseqs2 LCA outputs, then joins per-sample data into tables for ranks such as kingdom, phylum, class, order, family, genus, and species.

This is a deliberate extra step. TaxonKit can provide more complete lineages with fewer gaps than relying only on the native MMseqs2 output. It requires a current NCBI taxonomy dump and a correctly set `TAXONKIT_DB`.

The scheduler wrapper runs the Snakemake component and joining scripts. For standalone use, consult the repository's `summarise_taxonomy/README.md`.

## 6. Functional annotation with BV-BRC Subsystems

UniProt/UniRef identifiers from MMseqs2 are joined to an SQLite mapping derived from the BV-BRC Subsystems framework. Subsystems provide a hierarchical functional vocabulary. The workflow produces per-sample functional records and summary tables, including raw and normalised counts and taxonomic-functional combinations where supported.

Normalisation makes samples with different mapped-read depths more comparable, but it does not remove every source of technical or compositional bias. Preserve raw counts, document the denominator, and select downstream statistical methods appropriate for compositional data.

## 7. Optional PHROG annotation

Some paired profiles include scripts to add PHROG phage-protein functions and count them. This requires the relevant annotation/mapping data and is not required for the core taxonomy or Subsystems results.

## 8. Assembly with MEGAHIT

MEGAHIT assembles cleaned, non-host reads with a succinct de Bruijn graph. Profiles may assemble each sample independently to avoid the extreme memory demand of complete cross-assembly; other scripts support larger combined assemblies where cluster resources permit.

Assembly can run in parallel with read-based taxonomy once cleaned reads are ready. Review contig counts, total assembled bases, N50 alongside other statistics, minimum contig lengths, and warnings. An assembly completing successfully does not guarantee biological quality.

## 9. Preparing contigs for VAMB

VAMB jointly uses contig composition and differential abundance. Preparation generally includes:

1. concatenating assemblies while retaining unique sample/contig identity;
2. building a minimap2 index for the combined contigs; and
3. mapping every sample's reads back to the contigs to create abundance inputs.

The concatenation job must finish before mapping. All mapping array tasks must finish before VAMB runs. Check that sample names, BAM outputs, contig identifiers, and expected counts align.

## 10. Binning with VAMB

VAMB uses a variational autoencoder to combine k-mer composition and coverage information in a latent representation, then clusters contigs into candidate MAGs. Some profiles request GPUs and load a cluster-specific runtime such as ROCm; others may use different installations.

GPU directives, account suffixes, module versions, and VAMB builds are especially cluster-specific. Confirm them with a minimal GPU job. VAMB may emit split and unsplit clusters; keep both until quality assessment determines which set is useful.

## 11. MAG quality with CheckM

CheckM estimates completeness and contamination for candidate bins. Profiles can run it separately for split and unsplit VAMB clusters. Its reference data and software environment must be installed and visible to compute nodes.

Quality labels are screening criteria rather than proof that a bin represents a complete organism. Examine contamination, strain heterogeneity, taxonomy, coverage, and downstream biological consistency.

## 12. Optional grouped VAMB

Paired profiles can bin subsets independently using a two-column tab-separated mapping of group to sample prefix. The prefix must identify zero or one entry in `R1_reads.txt`; ambiguous substrings can silently select unintended samples.

Grouped processing repeats concatenation, mapping, and VAMB for each defined group. Use it when coverage patterns are expected to be more informative within biological or technical strata, and retain the mapping file as analysis provenance.

## 13. Read-fate and Sankey summaries

Where supplied, read-fate and Sankey scripts combine counts from quality control, host removal, annotation, and assembly-oriented stages. These summaries are valuable pipeline sanity checks: large unexplained losses often reveal naming, input, threshold, or failed-stage problems.
