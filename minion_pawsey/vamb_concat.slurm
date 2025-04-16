#!/bin/bash
#SBATCH --job-name=VambConcat
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH -o slurm_output/vamb_slurm/vamb_concat-%j.out
#SBATCH -e slurm_output/vamb_slurm/vamb_concat-%j.err

set -euo pipefail
if [[ ! -n "$ATAVIDE_CONDA" ]]; then
    echo "Please include the ATAVIDE_CONDA location on your slurm submission: --export=ATAVIDE_CONDA=$ATAVIDE_CONDA" >&2
    exit 2;
fi
eval "$(conda shell.bash hook)"
conda activate $ATAVIDE_CONDA


mkdir -p vamb
python ~/atavide_lite/bin/vamb_concatenate.py -m 300 vamb/contigs.fna.gz megahit/*/final.contigs.fa

