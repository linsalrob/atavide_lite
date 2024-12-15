#!/bin/bash
#SBATCH --job-name=fastq2fasta
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -o slurm_output/faing-%j.out
#SBATCH -e slurm_output/faing-%j.err

set -euo pipefail


if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh


# this requires ~/atavide_lite/bin/fastq2fasta.c be compiled with make, and it has built in reading gzip support

extract_R1 () {
	FQ=$1;
	FB=`basename $FQ .gz`;
	FA=${FB/fastq/fasta}
	~/atavide_lite/bin/fastq2fasta -n 1 $FQ - | gzip -c > fasta/$FA.gz;
}

export -f extract_R1
mkdir -p fasta

find -L $HOSTREMOVED -type f -name *fastq.gz | parallel -j 16 extract_R1

