#!/bin/bash
#SBATCH --job-name=ReadFate
#SBATCH --time=0-10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH -o slurm_output/read_fate-%j.out
#SBATCH -e slurm_output/read_fate-%j.err

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate atavide_lite

mkdir -p Unknowns

export PYTHONPATH=$PYTHONPATH:~/atavide_lite
python ~/atavide_lite/bin/read_fate.py -r reads.txt -u Unknowns/unknowns.fastq -v
