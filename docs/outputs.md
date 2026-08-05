# Outputs and interpretation

Exact filenames vary by profile and dataset configuration. This page describes the stable output concepts; the selected profile README and scripts define exact paths.

## Scheduler logs

`slurm_output/` and equivalent PBS locations contain standard output and error logs. Array filenames normally include both the parent job ID and array task ID. Preserve these logs at least until the analysis is validated because they record tool commands, versions, timings, warnings, and failure locations.

A missing downstream output should first be investigated through scheduler state and the corresponding `.err` file, not by resubmitting the entire workflow.

## Quality-controlled reads

`fastq_fastp/` commonly contains cleaned FASTQ files, while `fastp_output/` contains JSON and HTML reports. Review:

- reads before and after filtering;
- adapter trimming;
- quality by cycle;
- length distributions;
- duplication and GC patterns; and
- consistency between mates and samples.

## Host and non-host reads

Directories named by `HOST` and `HOSTREMOVED` contain the separated read sets. The non-host directory is the normal input to subsequent stages. Record mapping rates per sample and investigate outliers before annotation or assembly.

Host-mapped reads may be sensitive. Apply the appropriate access, retention, and deletion policy even if the pipeline does not use them downstream.

## MMseqs2 results

The `mmseqs/` hierarchy generally contains one directory per sample with `easy-taxonomy` results, including LCA and report files. Treat these as intermediate evidence supporting summary tables. Preserve the command, database release, and relevant MMseqs2 version/parameters.

## Taxonomy tables

Taxonomy processing creates per-record lineage files and joined tables at standard ranks. Typical summary outputs include kingdom, phylum, class, order, family, genus, and species tables, often compressed TSV.

Rows or columns depend on the joining script and workflow version. Before statistical analysis, confirm orientation, identifier uniqueness, unclassified handling, taxonomic naming conventions, and whether values are raw counts or another measure.

## Functional and Subsystems tables

Functional outputs join UniRef identifiers, protein descriptions, taxonomic information, and hierarchical Subsystems levels. Count summaries may include raw and normalised abundance.

A protein can participate in multiple functional assignments. Some utilities distribute a read's contribution across mappings using weighted counts. Read the output header and utility documentation rather than assuming every count is an integer.

## Assemblies

MEGAHIT creates one output hierarchy per sample in the relevant assembly directory. Important products include final contigs and, where generated, graph formats. MEGAHIT working directories can be large; copy final products to durable storage before cleaning temporary data.

Interpret assembly statistics together. A higher N50 alone does not establish a better metagenome assembly, and results depend strongly on depth, complexity, read quality, and host removal.

## VAMB products

VAMB-related directories include concatenated contigs, mapping products, abundance data, cluster assignments, and candidate bin FASTAs. Profiles may retain `clusters_split` and `clusters_unsplit`, group-specific directories, or legacy `bins/` paths.

Contig names are the link between assemblies, mappings, clusters, and bins. Do not rename them after concatenation unless all dependent records are changed consistently.

## CheckM results

CheckM directories contain bin-level completeness and contamination estimates. Keep the exact bin set and CheckM database/version associated with every report. Comparing reports from split and unsplit VAMB outputs can guide which binning strategy to carry forward.

## Sankey and read-fate outputs

These consolidate how many reads remain or receive assignments after key stages. Use them to detect:

- unexpected quality-control loss;
- unusually high host mapping;
- samples with few MMseqs2 assignments;
- mate-count inconsistencies; and
- stages that produced no counted output.

## Reproducibility record

For each completed analysis, archive or record:

- repository commit ID and selected profile;
- `DEFINITIONS.sh`, read list, and grouping files;
- scheduler submission commands and relevant logs;
- software environment exports and module versions;
- scheduler resources actually used;
- host and annotation database identities/checksums;
- tool parameters changed from profile defaults; and
- final output checksums and storage locations.
