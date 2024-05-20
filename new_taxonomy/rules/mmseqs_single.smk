###################################################
#                                                 #
# Run mmseqs easy taxonomy on a directory         #
# against the UniRef50 database                   #
#                                                 #
# This version works with minion (UNPAIRED)       #
# fasta files                                     #
#                                                 #
#                                                 #
#                                                 #
###################################################


import os
import sys

DB = "UniRef50"
DB_SOURCE = "/scratch/user/edwa0468/Databases/mmseqs/"
SCRATCH = os.environ['BGFS']
MMSEQS = "mmseqs"

if not os.path.exists(SCRATCH):
    print(f"ERROR: Scratch location {SCRATCH} does not exist", 
          file=sys.stderr)
    sys.exit(1)

if not os.path.exists(DB_SOURCE):
    print(f"ERROR: {DB_SOURCE} does not exist. Please set the location of your databases",
          file=sys.stderr)
    sys.exit(1)


FASTADIR = "fasta"

SAMPLES, = glob_wildcards(os.path.join(FASTADIR, "{sample}.fasta.gz"))

rule all:
    input:
        expand(os.path.join(MMSEQS, "{sample}", "{sample}_tophit_report.gz"), sample=SAMPLES)

rule run_mmseqs:
    input:
        fa = os.path.join(FASTADIR, "{sample}.fasta.gz"),
        db = os.path.join(DB_SOURCE, DB, DB)
    output:
        thr = os.path.join(MMSEQS, "{sample}", "{sample}_tophit_report.gz"),
        lca = os.path.join(MMSEQS, "{sample}", "{sample}_lca.tsv.gz"),
        rep = os.path.join(MMSEQS, "{sample}", "{sample}_report.gz"),
        aln = os.path.join(MMSEQS, "{sample}", "{sample}_tophit_aln.gz")
    params:
        scratch = SCRATCH,
        tmpout = os.path.join(SCRATCH, "{sample}"),
        finalout = os.path.join(MMSEQS, "{sample}")
    conda: "envs/mmseqs.yaml"
    resources:
        cores = 32,
        mem_mb = 200000
    shell:
        """
        mmseqs easy-taxonomy {input.fa} {input.db} {params.tmpout} $(mktemp -d -p {params.scratch})  --threads {resources.cores} &&
        rsync -avz {params.tmpout} 
        """

