#!/bin/bash
#SBATCH --job-name=vamb_split
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH -o slurm_output/vamb_slurm/%x-%j.out
#SBATCH -e slurm_output/vamb_slurm/%x-%j.err

set -euo pipefail

VENV_DIR="/scratch/pawsey1018/edwa0468/software/pip/vamb"
if [[ ! -e "$VENV_DIR" ]]; then 
	echo "Please run vamb_install.slurm to create the vamb virtual environment" >&2
	exit 1;
fi

module load rocm/6.2.4
source "$VENV_DIR/bin/activate"


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

