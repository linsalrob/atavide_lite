#!/bin/bash
#SBATCH --job-name=mmTax
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --partition=short
#SBATCH -o slurm_output/mmtax-%j.out
#SBATCH -e slurm_output/mmtax-%j.err


set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate atavide_lite

snakemake --profile slurm -s ~/GitHubs/atavide_lite/new_taxonomy/taxonomy.smk
python ~/GitHubs/atavide_lite/new_taxonomy/scripts/join_taxonomies.py -t taxonomy -o taxonomy_summary/
