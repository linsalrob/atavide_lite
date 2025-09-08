#!/bin/bash
#SBATCH --job-name=vamb_split
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH -o slurm_output/vamb_slurm/%x-%j.out
#SBATCH -e slurm_output/vamb_slurm/%x-%j.err

set -euo pipefail

if [[ ! -n "$ATAVIDE_CONDA" ]]; then
    echo "Please define the ATAVIDE_CONDA environment variable with the location of your ATAVIDE_LITE conda installation" >&2;
    exit 2;
fi

eval "$(conda shell.bash hook)"
conda activate $ATAVIDE_CONDA


if [[ $# -eq 0 ]]; then 
	echo "$0 <vamb contigs> <clusters file> <min bin size> <output directory>"
	exit;
fi

CONTIGS=$1
CLUSTERS=$2
MINBIN=$3
ODIR=$4


echo "Running: ~/atavide_lite/bin/vamb_create_fasta.py $CONTIGS $CLUSTERS $MINBIN $ODIR" >&2

python ~/atavide_lite/bin/vamb_create_fasta.py $CONTIGS $CLUSTERS $MINBIN $ODIR

