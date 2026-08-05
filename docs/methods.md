# Scientific methods

This page describes the scientific intent of the shipped workflow. Exact commands and defaults remain visible in each cluster profile and should be recorded for every analysis.

## Quality control

Read-based analyses begin with adapter trimming and filtering using fastp. Profile defaults are tailored separately for paired short reads and single/long reads. Paired profiles commonly retain reads of at least 100 bp with no more than one uncalled base, while long-read examples may use different thresholds. Users should select thresholds appropriate to the library and preserve fastp reports.

## Host-read removal

Quality-controlled reads are aligned to a user-selected host reference with minimap2. Samtools separates mapped and unmapped primary reads that pass platform/vendor checks and are not marked as PCR or optical duplicates. Paired profiles preserve R1/R2 identity through additional flag selection. Non-host reads continue through metagenomics analysis; host-mapped reads are stored separately but are not currently used downstream.

The host reference can contain any genome or multi-FASTA sequence collection that the analysis should exclude. Its accession, version, and checksum are part of the method.

## Read-based taxonomic profiling

Cleaned reads are converted to FASTA using the repository's compiled converter and analysed with MMseqs2 `easy-taxonomy`. Current project workflows favour UniRef100 for improved hit coverage, although UniRef50 variants remain available for comparison and lower resource demand.

MMseqs2 uses k-mer prefiltering, ungapped filtering, vectorised alignment, and a lowest-common-ancestor procedure to assign taxonomy. Results depend on the MMseqs2 and UniRef releases and any changed parameters.

Taxonomic lineages are expanded using TaxonKit through PyTaxonKit and merged into per-rank tables. This additional step was chosen because it can recover more complete standard lineages than the native MMseqs2 representation alone.

## Functional profiling

UniRef identifiers from MMseqs2 results are joined to protein functions and the hierarchical BV-BRC Subsystems framework through a project SQLite mapping. Mapped reads are counted for each function/Subsystem. Outputs can include raw counts, counts normalised by mapped-read totals, and combined taxonomic-functional records. Where one protein maps to multiple functions, weighted counts may divide its contribution among assignments.

## Assembly

Cleaned non-host reads are assembled with MEGAHIT, which uses a succinct de Bruijn graph across multiple k-mer sizes and manages memory conservatively. Profiles commonly assemble samples independently to avoid the memory growth of complete cross-assembly, while specialised scripts support larger combined assemblies where appropriate resources are available.

## MAG recovery

Assemblies are concatenated with unique contig identifiers. Minimap2 maps reads from each sample back to the combined contigs to estimate differential abundance. VAMB combines abundance with tetranucleotide composition in a variational-autoencoder representation and clusters contigs into candidate MAGs. Depending on the VAMB version/profile, split and unsplit cluster sets may be produced.

CheckM estimates completeness and contamination for candidate bins. These estimates support quality screening but should be interpreted with taxonomy, coverage, strain variation, and the intended downstream use.

## Optional analyses

Profiles may provide PHROG functional annotation for phage proteins, grouped VAMB runs for user-defined sample strata, read-fate/Sankey summaries, 16S screening, and selective removal of taxonomic lineages from downstream summaries.

## Reproducible reporting

A methods description should identify:

- the atavide lite commit and profile;
- sequencing type and sample naming model;
- quality-control parameters and adapters;
- host reference and filtering decision;
- MMseqs2, UniRef, and taxonomy releases;
- Subsystems mapping release;
- MEGAHIT and VAMB parameters;
- split/unsplit/grouped binning choice;
- CheckM version/database; and
- any changes to supplied scheduler resources or commands.

Use the repository's `citation.cff` to cite atavide lite and `references.bib` to identify the primary publications for integrated tools.
