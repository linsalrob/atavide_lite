################################################################
#                                                              #
# ORF-based read taxonomy.                                     #
#                                                              #
# NOTE: This is designed to be used with minion reads          #
# therefore, we don't have R1/R2, just fastq files             #
#                                                              #
# (c) Rob Edwards                                              #
#                                                              #
#                                                              #
#                                                              #
################################################################

import os
import sys

if not os.environ["TAXONKIT_DB"]:
    print("Please define the TAXONKIT_DB database location (even if you set it to ~/.taxonkit)", file=sys.stderr)
    sys.exit(1)


READDIR = 'fastq'
ORFDIR  = 'orfs'
FASTADIR = ORFDIR
MMSEQS = 'mmseqs'
TAX = 'taxonomy'
SUBSYSTEMS = "subsystems"

# A Snakemake regular expression matching the forward mate FASTQ files.
# the comma after SAMPLES is important!
FQSAMPLES, = glob_wildcards(os.path.join(READDIR, '{sample}.fastq.gz'))
if len(FQSAMPLES) == 0:
    sys.stderr.write(f"We did not find any fastq files in {SAMPLES}. Is this the right read dir?\n")
    sys.exit(0)


#localrules: all_contig

# include the required rules
include: "rules/get_orfs.smk"
include: "rules/mmseqs_single.smk"
include: "rules/taxonomy.smk"
include: "rules/summarise.smk"
include: "rules/add_functions.smk"
include: "rules/summarise_subsystem_counts.smk"


# make some output directories
os.makedirs(ORFDIR, exist_ok=True)
os.makedirs(MMSEQS, exist_ok=True)

rule all_contig:
    input:
        expand(os.path.join(MMSEQS, "{sample}", "{sample}.contigs.tsv.gz"), sample=FQSAMPLES),
        expand(os.path.join(ORFDIR, '{sample}.fasta.gz'), sample=FQSAMPLES),
        expand(os.path.join(MMSEQS, "{sample}", "{sample}_lca.tsv.gz"), sample=FQSAMPLES),
        expand(os.path.join(MMSEQS, "{sample}", "{sample}_lca_taxonomy.tsv.gz"), sample=FQSAMPLES),
        expand(os.path.join(MMSEQS, "{sample}", "{sample}_tophit_report_subsystems.gz"), sample=FQSAMPLES),
        os.path.join(TAX, "kingdom.tsv.gz"), os.path.join(TAX, "phylum.tsv.gz"), 
        os.path.join(TAX, "class.tsv.gz"), os.path.join(TAX, "order.tsv.gz"), 
        os.path.join(TAX, "family.tsv.gz"), os.path.join(TAX, "genus.tsv.gz"), 
        os.path.join(TAX, "species.tsv.gz"), 
        SUBSYSTEMS

rule concatentate_taxonomies:
    input:
        taxfiles = expand(os.path.join(MMSEQS, "{sample}", "{sample}.contigs.tsv.gz"), sample=FQSAMPLES)
    output:
        k = os.path.join(TAX, "kingdom.tsv.gz"),
        p = os.path.join(TAX, "phylum.tsv.gz"),
        c = os.path.join(TAX, "class.tsv.gz"),
        o = os.path.join(TAX, "order.tsv.gz"),
        f = os.path.join(TAX, "family.tsv.gz"),
        g = os.path.join(TAX, "genus.tsv.gz"),
        s = os.path.join(TAX, "species.tsv.gz"),
        outdir = directory(TAX)
    script:
        "scripts/combine_taxonomies.py"

