#!/bin/bash
#SBATCH --job-name=VambConcat
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH -o slurm_output/vamb_slurm/vamb_concat-%j.out
#SBATCH -e slurm_output/vamb_slurm/vamb_concat-%j.err

set -euo pipefail

#VENV_DIR="/scratch/$PAWSEY_PROJECT/$USER/software/pip/vamb"
VENV_DIR="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/vamb"
if [[ ! -x "$VENV_DIR/bin/vamb" ]]; then
	echo "Please run vamb_install.slurm to create the vamb virtual environment" >&2
	exit 1;
fi

# source "$VENV_DIR/bin/activate"
eval "$(conda shell.bash hook)"
conda activate $VENV_DIR


mkdir -p vamb_crass
python ~/GitHubs/atavide_lite/bin/vamb_concatenate.py -m 300 vamb_crass/contigs.fna.gz megahit_all_reads/final.contigs.fa

