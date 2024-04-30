# a new way to calculate the taxonomy files, using taxontoolkit

import os
import sys

if not os.environ["TAXONKIT_DB"]:
    print("Please define the TAXONKIT_DB database location (even if you set it to ~/.taxonkit)", file=sys.stderr)
    sys.exit(1)


MMSEQS = "mmseqs"
TEST="SAGCFN_22_01149_S3"
SAMPLES, LCAS, = glob_wildcards(os.path.join(MMSEQS, '{sample}', '{lca_samp}_lca.tsv.gz'))


rule all:
    input:
        expand(os.path.join(MMSEQS, "{sample}", "{sample}_lca_taxonomy.tsv.gz"), sample=SAMPLES)

rule list_taxonomy:
    input:
        lca = os.path.join(MMSEQS,  "{sample}", "{sample}_lca.tsv.gz")
    output:
        tk = os.path.join(MMSEQS, "{sample}", "{sample}_taxonomy.tsv.gz")
    conda: "envs/taxonkit.yaml"
    resources:
        cores=4,
        mem_mb=16000
    shell:
        """
        zcat {input.lca} | awk '!s[$2]++ {{print $2}}' | taxonkit lineage | taxonkit reformat -P | awk -F"\t" -v OFS="\t" '{{print $1,$3}}' | gzip -c > {output.tk} 
        """


rule add_taxonomy:
    input:
        tk = os.path.join(MMSEQS, "{sample}", "{sample}_taxonomy.tsv.gz"),
        lca = os.path.join(MMSEQS,  "{sample}", "{sample}_lca.tsv.gz")
    output:
        lcatax = os.path.join(MMSEQS,  "{sample}", "{sample}_lca_taxonomy.tsv.gz")
    resources:
        mem_mb=8000
    script: "scripts/merge_taxonomy.py"


