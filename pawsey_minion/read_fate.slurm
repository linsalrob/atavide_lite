#!/bin/bash
#SBATCH --job-name=ReadFate
#SBATCH --time=0-10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH -o slurm_output/read_fate-%j.out
#SBATCH -e slurm_output/read_fate-%j.err

set -euo pipefail
if [[ ! -n "$ATAVIDE_CONDA" ]]; then
    echo "Please include the ATAVIDE_CONDA location on your slurm submission: --export=ATAVIDE_CONDA=$ATAVIDE_CONDA" >&2
    exit 2;
fi
eval "$(conda shell.bash hook)"
conda activate $ATAVIDE_CONDA

mkdir -p Unknowns

export PYTHONPATH=$PYTHONPATH:~/atavide_lite
python ~/atavide_lite/bin/read_fate.py -r reads.txt -u Unknowns/unknowns.fastq -v
