################################################################
#                                                              #
# Get the ORFs from a fastq file and put them in a fasta file  #
#                                                              #
#                                                              #
# (c) Rob Edwards                                              #
#                                                              #
#                                                              #
#                                                              #
################################################################

import os
import sys


rule get_orfs:
    """
    Get the open reading frames. Note that this is hard coded at the momnet
    because bioconda is broken, and I can't add anything!
    """

    input:
        fq = os.path.join(READDIR, '{sample}.fastq.gz')
    output:
        faa = os.path.join(ORFDIR, '{sample}.fasta')
    params:
        odir = ORFDIR,
        minorflen = 30
    resources:
        cores = 16
    conda: "../envs/get_orfs.yaml"
    shell:
        """
        get_orfs -f {input.fq}  -l {params.minorflen} -j {resources.cores} > {output.faa}
        """

rule compress_fastq:
    input:
        faa = os.path.join(ORFDIR, '{sample}.fasta')
    output:
        faa = os.path.join(ORFDIR, '{sample}.fasta.gz')
    resources:
        cores = 8
    shell:
        """
        pigz {input.faa}
        """

