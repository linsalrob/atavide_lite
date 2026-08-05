# Components and utilities

The cluster profiles coordinate tools in the repository root and several component directories. This page helps users locate the implementation behind a stage.

## `bin/`

`bin/` contains compiled C programs, Python utilities, Perl, and shell helpers. Run `make all` there before using compiled tools.

Important utilities include:

- `fastq2fasta`: compiled FASTQ-to-FASTA conversion used before MMseqs2;
- `fastg2gfa`: converts MEGAHIT graph output;
- `vamb_concatenate.py`: prepares uniquely named concatenated contigs for VAMB;
- `vamb_create_fasta.py` and related scripts: materialise cluster/bin FASTAs;
- `count_subsystems.py` and selective variants: build functional count summaries;
- `taxonomy_selected.py`: filter downstream taxonomic summaries by a selected lineage;
- `read_fate.py` and `sankey_plot.py`: summarise pipeline disposition of reads;
- MMseqs2 report conversion and taxonomy/function joining utilities; and
- `sqlite/` scripts for building the MMseqs2/Subsystems mapping database.

Most Python utilities provide command-line help:

```bash
python bin/count_subsystems.py --help
```

Run help before reusing a utility outside its profile, because input columns and output assumptions are specialised.

## `summarise_taxonomy/`

This Snakemake component enriches MMseqs2 LCA output with TaxonKit lineages, merges sample results, and builds rank-specific summaries. It requires NCBI taxonomy dump files and `TAXONKIT_DB`.

The cluster wrapper usually runs this component. For standalone execution, the component README documents its expected `mmseqs/<sample>/` layout and joining command.

## `minion_taxonomy_and_function/`

This component supports a protein/ORF-oriented long-read taxonomy and function workflow. It assumes basecalling, long-read quality control, and host removal have already occurred, then predicts ORFs, runs amino-acid MMseqs2 taxonomy, completes lineages, and adds BV-BRC Subsystems functions.

It is more specialised than the main `pawsey_minion`/`deepthought_minion` script sequence. Read its README and Snakefiles before deciding which long-read route fits the analysis.

## `adapters/`

Adapter FASTA files are provided for paired short-read and Oxford Nanopore library preparation examples. Confirm the library kit and adapter set match the data. Kit names and chemistry change over time.

## `pawsey_lib/`

Pawsey helper scripts check or recreate temporary Conda environments expected by Pawsey profiles. They encode Setonix paths and should not be reused unchanged on another cluster.

## `bash/`

The single Bash workflow is an educational, linear representation of the major commands. It is helpful when understanding or porting the pipeline, but running stages independently preserves the project's failure isolation and scheduler-resource benefits.

## Custom taxon exclusion

After taxonomy and Subsystems processing, selected lineages can be removed to test whether a signal depends on them. For example:

```bash
python bin/count_subsystems_selective.py \
  --directory mmseqs \
  --subsystems subsystems_no_pseudo \
  --name no_pseudo \
  --taxa f__Pseudomonadaceae \
  --verbose
```

Then create a filtered taxonomy hierarchy and summaries with `taxonomy_selected.py` and `join_taxonomies.py`. Check each script's current `--help`; option names may differ between versions. Preserve both original and filtered tables and document the selected lineage.
