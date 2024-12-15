#!/bin/bash
#SBATCH --job-name=VambConcat
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH -o slurm_output/vamb_slurm/vamb_concat-%j.out
#SBATCH -e slurm_output/vamb_slurm/vamb_concat-%j.err

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate /scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite_vamb

mkdir -p vamb
python ~/atavide_lite/bin/vamb_concatenate.py -m 300 vamb/contigs.fna.gz megahit/*/final.contigs.fa

