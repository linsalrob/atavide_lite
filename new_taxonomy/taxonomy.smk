# a new way to calculate the taxonomy files, using taxontoolkit

import os
import sys

if not os.environ["TAXONKIT_DB"]:
    print("Please define the TAXONKIT_DB database location (even if you set it to ~/.taxonkit)", file=sys.stderr)
    sys.exit(1)


MMSEQS = "mmseqs"
TEST="SAGCFN_22_01149_S3"
SAMPLES, LCAS, = glob_wildcards(os.path.join(MMSEQS, '{sample}', '{lca_samp}_lca.tsv.gz'))
TAXDIR="taxonomy"
TAXOUTPUTDIR="taxonomy_summary"
os.makedirs(TAXDIR, exist_ok=True)
os.makedirs(TAXOUTPUTDIR, exist_ok=True)
TAXONOMIES = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

# we have a checkpoint for the summarised files, so this generates a list of those files
def summarised_files(wildcards):
    ck_output = checkpoints.summarise.get(**wildcards).output.directory
    print("Reading directory {ck_output}", file-sys.stderr)
    SMP,TXS, = glob_wildcards(os.path.join(ck_output, "{smps}", "{tx}.tsv.gz"))
    return expand(os.path.join(ck_output, "{smp}", f"{tax}.tsv.gz"), smp=SMP, tax=TXS)

# files = os.path.join(TAXDIR, "{sample}", "{tax}.tsv.gz"),
# # Use a function to dynamically generate the input files based on the checkpoint
#         lambda wildcards: expand("data/{sample}_{tax}.txt", sample=checkpoint_output("collect_files")["samples"], tax=wildcards.tax)


rule all:
    input:
        expand(
            [os.path.join(MMSEQS, "{sample}", "{sample}_lca_taxonomy.tsv.gz"),
            os.path.join(TAXDIR, "{sample}")], 
            sample=SAMPLES),
        # expand(os.path.join(TAXDIR, "{tax}.tsv.gz"), tax=TAXONOMIES)


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



checkpoint summarise:
    input:
        lcatax = os.path.join(MMSEQS, "{sample}", "{sample}_lca_taxonomy.tsv.gz")
    output:
        directory = directory(os.path.join(TAXDIR, "{sample}"))
    script:
        "scripts/summarise_taxonomy.py"

"""
rule combine:
    input:
        files = summarised_files
    output:
        raw_output = os.path.join(TAXDIR, "{tax}.raw.tsv.gz"),
        norm_output = os.path.join(TAXDIR, "{tax}.norm.tsv.gz")
    script:
        "scripts/join.py"
"""

rule combine_all:
    input:
        summarised_files
    output:
        os.path.join(TAXOUTPUTDIR, "all_taxonomies.tsv")
    shell:
        """
        cat {input} > {output}
        """
