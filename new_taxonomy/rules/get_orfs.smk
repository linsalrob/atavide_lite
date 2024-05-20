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


READDIR = 'fastq'
ORFDIR  = 'orfs'

# A Snakemake regular expression matching the forward mate FASTQ files.
# the comma after SAMPLES is important!
SAMPLES, = glob_wildcards(os.path.join(READDIR, '{sample}.fastq.gz'))
if len(SAMPLES) == 0:
    sys.stderr.write(f"We did not find any fastq files in {SAMPLES}. Is this the right read dir?\n")
    sys.exit(0)

os.makedirs(ORFDIR, exist_ok=True)

rule all:
    input:
        expand(os.path.join(ORFDIR, '{sample}.faa.gz'), sample=SAMPLES)

rule get_orfs:
    """
    Get the open reading frames. Note that this is hard coded at the momnet
    because bioconda is broken, and I can't add anything!
    """

	input:
        fq = os.path.join(READDIR, '{sample}.fastq.gz')
    output:
        faa = os.path.join(ORFDIR, '{sample}.fasta.gz')
    params:
        odir = ORFDIR,
        minorflen = 30
    shell:
        """
        ~/GitHubs/get_orfs/build/get_orfs -f {input.fq}  -l {params.minorflen} > {output.faa}
        """

