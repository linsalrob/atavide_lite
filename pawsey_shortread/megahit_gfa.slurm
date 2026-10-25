#!/bin/bash
#SBATCH --job-name=megahit_gfa
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1800M
#SBATCH -o slurm_output/megahit_slurm/megahit_gfa-%j.out
#SBATCH -e slurm_output/megahit_slurm/megahit_gfa-%j.err

# NOTE: Setonix enforces MaxMemPerCPU (work=1840M, highmem=7900M). A bare --mem larger
# than cpus-per-task * MaxMemPerCPU does not fail - Slurm silently raises the CPU count
# to cover it. Request memory with --mem-per-cpu so the two can never disagree.


set -euo pipefail
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# test we have atavide_lite conda installed, and then activate it
eval "$(conda shell.bash hook)"
"$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh"
ATAVIDE_CONDA="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite"
conda activate $ATAVIDE_CONDA

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file. This script optionally reads ASSDIR and KMER from it" >&2
	exit 2;
fi

source DEFINITIONS.sh

# Allow DEFINITIONS.sh (sourced above) or the environment to override these; the defaults
# preserve the previous hardcoded behaviour.
ASSDIR=${ASSDIR:-megahit_all_reads}
KMER=${KMER:-141}

echo "Converting fa to fastg"
megahit_toolkit contig2fastg $KMER $ASSDIR/intermediate_contigs/k$KMER.contigs.fa> $ASSDIR/k$KMER.fastg
echo "Converting fastg to gfa"
agtools fastg2gfa -g $ASSDIR/k$KMER.fastg -k $KMER -o $ASSDIR/k$KMER.gfa
echo "Fin"

