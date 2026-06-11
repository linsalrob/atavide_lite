#!/bin/bash
#SBATCH --job-name=fastq2fasta
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH -o slurm_output/faing-%j.out
#SBATCH -e slurm_output/faing-%j.err


set -euo pipefail
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# test we have atavide_lite conda installed, and then activate it
eval "$(conda shell.bash hook)"
"$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh"
ATAVIDE_CONDA="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite"
conda activate $ATAVIDE_CONDA

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

extract_R2 () {
	FQ=$1;
	FB=`basename $FQ .gz`;
	FA=${FB/fastq/fasta}
	~/atavide_lite/bin/fastq2fasta -n 2 $FQ - | gzip -c > fasta/$FA.gz;
}

export -f extract_R1
export -f extract_R2
mkdir -p fasta

find -L $HOSTREMOVED -type f -name \*_R1\*fastq.gz | parallel -j 16 extract_R1
find -L $HOSTREMOVED -type f -name \*_R2\*fastq.gz | parallel -j 16 extract_R2

