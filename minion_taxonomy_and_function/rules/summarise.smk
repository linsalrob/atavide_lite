


rule summarise:
    input:
        lcatax = os.path.join(MMSEQS, "{sample}", "{sample}_lca_taxonomy.tsv.gz")
    output:
        consum = os.path.join(MMSEQS, "{sample}", "{sample}.contigs.tsv.gz")
    script:
        "../scripts/summarise_contig_taxonomy.py"
