#!/bin/bash
#SBATCH --job-name=VambConcat
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH -o slurm_output/vamb_slurm/vamb_concat-%j.out
#SBATCH -e slurm_output/vamb_slurm/vamb_concat-%j.err

set -euo pipefail

if [[ ! -n "$ATAVIDE_CONDA_VAMB" ]]; then
    echo "Please define the ATAVIDE_CONDA_VAMB environment variable with the location of your ATAVIDE_LITE_VAMB conda installation" >&2;
    echo "Note: this is based on the atavide_vamb.yaml and is different from the regular atavide conda" >&2;
    exit 2;
fi

eval "$(conda shell.bash hook)"
conda activate $ATAVIDE_CONDA_VAMB

mkdir -p vamb
python ~/atavide_lite/bin/vamb_concatenate.py -m 300 vamb/contigs.fna.gz megahit/*/output/final.contigs.fa

