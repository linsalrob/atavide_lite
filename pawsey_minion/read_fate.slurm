#!/bin/bash
#SBATCH --job-name=ReadFate
#SBATCH --time=0-10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH -o slurm_output/read_fate-%j.out
#SBATCH -e slurm_output/read_fate-%j.err

set -euo pipefail
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# test we have atavide_lite conda installed, and then activate it
eval "$(conda shell.bash hook)"
"$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh"
ATAVIDE_CONDA="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite"
conda activate $ATAVIDE_CONDA

mkdir -p Unknowns

export PYTHONPATH=$PYTHONPATH:~/atavide_lite
python ~/atavide_lite/bin/read_fate.py -r reads.txt -u Unknowns/unknowns.fastq -v
