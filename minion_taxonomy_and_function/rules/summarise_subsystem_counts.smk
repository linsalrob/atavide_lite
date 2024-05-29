
rule summarise_ss_counts:
    input:
        ss = expand(os.path.join(MMSEQS, "{sample}", "{sample}_tophit_report_subsystems.gz"), sample=FQSAMPLES)
    output:
        odir = directory(SUBSYSTEMS)
    script: "../scripts/count_subsystems.py"

