#!/bin/bash
#SBATCH --job-name=VambConcat
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=short
#SBATCH -o vamb_concat-%j.out
#SBATCH -e vamb_concat-%j.err

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate vamb

mkdir -p vamb
python /home/edwa0468/GitHubs/atavide_lite/bin/vamb_concatenate.py vamb/contigs.fna.gz megahit/*/output/final.contigs.fa

