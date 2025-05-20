#!/bin/bash
#SBATCH --account=pawsey1018
#SBATCH --job-name=sankey_plot
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64GB
#SBATCH -o sankey_plot-%j.out
#SBATCH -e sankey_plot-%j.err

set -euo pipefail

export PYTHONPATH=$PYTHONPATH:~/GitHubs/atavide_lite
python ~/GitHubs/atavide_lite/bin/sankey_plot.py -r R1_reads.txt  --paired -v

