#!/bin/bash
#SBATCH --job-name=fastq2fasta
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH -o slurm_output/megahit_gfa-%j.out
#SBATCH -e slurm_output/megahit_gfa-%j.err


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

ASSDIR=megahit_all_reads
KMER=141

echo "Converting fa to fastg"
megahit_core contig2fastg $KMER $ASSDIR/intermediate_contigs/k$KMER.contigs.fa> $ASSDIR/k$KMER.fastg
echo "Converting fastg to gfa"
agtools fastg2gfa -g $ASSDIR/k$KMER.fastg -k $KMER -o $ASSDIR/k$KMER.gfa
echo "Fin"

