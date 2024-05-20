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
MMSEQS = 'mmseqs'

# A Snakemake regular expression matching the forward mate FASTQ files.
# the comma after SAMPLES is important!
SAMPLES, = glob_wildcards(os.path.join(READDIR, '{sample}.fastq.gz'))
if len(SAMPLES) == 0:
    sys.stderr.write(f"We did not find any fastq files in {SAMPLES}. Is this the right read dir?\n")
    sys.exit(0)

# include the required rules
include: "rules/get_orfs.smk"
include: "rules/mmseqs_single.smk"
include: "taxonomy.smk"
include: "rules/summarise.smk"

# make some output directories
os.makedirs(ORFDIR, exist_ok=True)
os.makedirs(MMSEQS, exist_ok=True)

rule all:
    input:
        expand(os.path.join(ORFDIR, '{sample}.fasta.gz'), sample=SAMPLES),
        expand(os.path.join(MMSEQS, "{sample}", "{sample}_lca.tsv.gz"), sample=SAMPLES),
        expand(os.path.join(MMSEQS, "{sample}", "{sample}_lca_taxonomy.tsv.gz"), sample=SAMPLES),
        expand(os.path.join(MMSEQS, "{sample}", "{sample}.contigs.tsv.gz"), sample=SAMPLES)



