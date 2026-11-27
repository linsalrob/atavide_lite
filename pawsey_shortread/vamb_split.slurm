#!/bin/bash
#SBATCH --job-name=vamb_split
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH -o slurm_output/vamb_slurm/%x-%j.out
#SBATCH -e slurm_output/vamb_slurm/%x-%j.err

set -euo pipefail

VENV_DIR="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/vamb"
if [[ ! -x "$VENV_DIR/bin/vamb" ]]; then
	echo "Please run vamb_install.slurm to create the vamb virtual environment" >&2
	exit 1;
fi

eval "$(conda shell.bash hook)"
conda activate "$VENV_DIR"


if [[ $# -eq 0 ]]; then 
	echo "$0 <vamb contigs> <clusters file> <min bin size> <output directory>"
	exit;
fi

CONTIGS=$1
CLUSTERS=$2
MINBIN=$3
ODIR=$4


echo "Running: ~/GitHubs/atavide_lite/bin/vamb_create_fasta.py $CONTIGS $CLUSTERS $MINBIN $ODIR" >&2

python ~/GitHubs/atavide_lite/bin/vamb_create_fasta.py $CONTIGS $CLUSTERS $MINBIN $ODIR

find $ODIR -type f -name "*.fna" -exec pigz {} \;
