#!/bin/bash
#SBATCH --job-name=mmTax
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH -o slurm_output/mmtax-%j.out
#SBATCH -e slurm_output/mmtax-%j.err


set -euo pipefail
if [[ ! -n "$ATAVIDE_CONDA" ]]; then
    echo "Please define the ATAVIDE_CONDA environment variable with the location of your ATAVIDE_LITE conda installation" >&2;
    exit 2;
fi

eval "$(conda shell.bash hook)"
conda activate $ATAVIDE_CONDA

snakemake --profile slurm -s ~/GitHubs/atavide_lite/new_taxonomy/taxonomy.smk
python ~/GitHubs/atavide_lite/new_taxonomy/scripts/join_taxonomies.py -t taxonomy -o taxonomy_summary/
