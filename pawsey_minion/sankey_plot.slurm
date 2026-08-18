#!/bin/bash
#SBATCH --job-name=sankey_plot
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64GB
#SBATCH -o slurm_output/sankey_plot-%j.out
#SBATCH -e slurm_output/sankey_plot-%j.err

set -euo pipefail

QC_DIRECTORY=${QC_DIR:-fastq_fastp}
QC_LABEL=${QC_LABEL:-fastp}

export PYTHONPATH=$PYTHONPATH:~/GitHubs/atavide_lite
python ~/GitHubs/atavide_lite/bin/sankey_plot.py -r reads.txt -q "$QC_DIRECTORY" --qc-label "$QC_LABEL" -t "$SLURM_CPUS_PER_TASK" -v
