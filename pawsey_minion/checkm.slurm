#!/bin/bash
#SBATCH --job-name=checkm
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128GB
#SBATCH -o slurm_output/checkm-%j.out
#SBATCH -e slurm_output/checkm-%j.err


set -euo pipefail

if [[ ! -n "$ATAVIDE_CONDA" ]]; then
    echo "Please include the ATAVIDE_CONDA location on your slurm submission: --export=ATAVIDE_CONDA=$ATAVIDE_CONDA" >&2
    exit 2;
fi
eval "$(conda shell.bash hook)"
conda activate $ATAVIDE_CONDA


if [[ -z $CHECKM_DATA_PATH ]]; then 
	echo "Please define CHECKM_DATA_PATH with your databases" >&2; 
	exit 1;
fi

if [[ $# -lt 2 ]]; then
	echo "$0 <source directory> <destination directory>" >&2
	exit 1
fi

# this %/ removes the trailing / if there is one
SOURCE=${1%/}
DEST=${2%/}
LOCALSOURCE=`basename $SOURCE`


echo "Running checkm on $LOCALSOURCE and output in $DEST"  >&2;
checkm lineage_wf -t 16 $LOCALSOURCE $DEST/checkm
checkm qa -t 16 $DEST/checkm/lineage.ms $DEST/checkm

