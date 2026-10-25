#!/bin/bash
#SBATCH --job-name=sankey_plot
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1800M
#SBATCH -o slurm_output/sankey_plot-%j.out
#SBATCH -e slurm_output/sankey_plot-%j.err

# NOTE: Setonix enforces MaxMemPerCPU (work=1840M, highmem=7900M). A bare --mem
# larger than cpus-per-task * MaxMemPerCPU does not fail - Slurm silently raises the
# CPU count to cover it, so the job queues for far more CPUs than it uses. Request
# memory with --mem-per-cpu so the two can never disagree.

set -euo pipefail

QC_DIRECTORY=${QC_DIR:-fastq_fastp}
QC_LABEL=${QC_LABEL:-fastp}

export PYTHONPATH=$PYTHONPATH:~/GitHubs/atavide_lite
python ~/GitHubs/atavide_lite/bin/sankey_plot.py -r reads.txt -q "$QC_DIRECTORY" --qc-label "$QC_LABEL" -t "$SLURM_CPUS_PER_TASK" -v
