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
DB_SOURCE = "/home/edwa0468/Databases/mmseqs"
MMSEQS = "mmseqs"
SCRATCH = os.environ["BGFS"]


if not os.path.exists(DB_SOURCE):
    print(f"ERROR: {DB_SOURCE} does not exist. Please set the location of your databases",
          file=sys.stderr)
    sys.exit(1)

if not os.path.exists(SCRATCH):
    os.makedirs(SCRATCH, exist_ok=True)


rule run_mmseqs:
    input:
        fa = os.path.join(FASTADIR, "{sample}.fasta.gz"),
        db = os.path.join(DB_SOURCE, DB, DB)
    output:
        thr = os.path.join(MMSEQS, "{sample}", "{sample}_tophit_report"),
        lca = os.path.join(MMSEQS, "{sample}", "{sample}_lca.tsv"),
        rep = os.path.join(MMSEQS, "{sample}", "{sample}_report"),
        aln = os.path.join(MMSEQS, "{sample}", "{sample}_tophit_aln")
    params:
        scratch = SCRATCH,
        tmpout = os.path.join(SCRATCH, "{sample}"),
        mmseqs = MMSEQS
    conda: "../envs/mmseqs.yaml"
    resources:
        cores = 32,
        mem_mb = 200000
    shell:
        """
        mkdir -p {params.tmpout} && 
        mmseqs easy-taxonomy {input.fa} {input.db} {params.tmpout}/{wildcards.sample} $(mktemp -d -p {params.scratch})  --threads {resources.cores} &&
        rsync -avz {params.tmpout} {params.mmseqs}
        """


rule compress_mmseqs:
    input:
        thr = os.path.join(MMSEQS, "{sample}", "{sample}_tophit_report"),
        lca = os.path.join(MMSEQS, "{sample}", "{sample}_lca.tsv"),
        rep = os.path.join(MMSEQS, "{sample}", "{sample}_report"),
        aln = os.path.join(MMSEQS, "{sample}", "{sample}_tophit_aln")
    output:
        thr = os.path.join(MMSEQS, "{sample}", "{sample}_tophit_report.gz"),
        lca = os.path.join(MMSEQS, "{sample}", "{sample}_lca.tsv.gz"),
        rep = os.path.join(MMSEQS, "{sample}", "{sample}_report.gz"),
        aln = os.path.join(MMSEQS, "{sample}", "{sample}_tophit_aln.gz")
    resources:
        cores = 8
    shell:
        """
        for F in {input}; do pigz $F; done
        """
