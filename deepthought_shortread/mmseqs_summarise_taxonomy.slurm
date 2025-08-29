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

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh

snakemake --profile slurm -s ~/GitHubs/atavide_lite/summarise_taxonomy/taxonomy.smk
python ~/GitHubs/atavide_lite/summarise_taxonomy/scripts/join_taxonomies.py -t taxonomy -o taxonomy_summary/ -n $SAMPLENAME
