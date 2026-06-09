#!/bin/bash
#SBATCH --job-name=checkm
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128GB
#SBATCH --partition=highmem
#SBATCH -o slurm_output/checkm-%j.out
#SBATCH -e slurm_output/checkm-%j.err

if [[ ! -n "$ATAVIDE_CONDA" ]]; then
    echo "Please define the ATAVIDE_CONDA environment variable with the location of your ATAVIDE_LITE conda installation" >&2;
    exit 2;
fi

if [[ -z $CHECKM_DATA_PATH ]]; then 
	CHECKM_DATA_PATH="$ATAVIDE_CONDA/checkm_data"
	if [[ ! -e $CHECKM_DATA_PATH ]]; then
		echo "Please install the checkm data (e.g. in $CHECKM_DATA_PATH) or set the variable CHECKM_DATA_PATH" >&2;
		exit 1;
	fi
	echo "Using $CHECKM_DATA_PATH for checkm data!"
fi


if [[ $# -lt 2 ]]; then
	echo "$0 <source directory> <destination directory>" >&2
	exit 1
fi

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate $ATAVIDE_CONDA

# this %/ removes the trailing / if there is one
SOURCE=${1%/}
DEST=${2%/}
LOCALSOURCE=`basename $SOURCE`


echo "Running checkm on $SOURCE and output in checkm"  >&2;
checkm lineage_wf -t $SLURM_CPUS_PER_TASK $SOURCE $DEST
checkm qa -t $SLURM_CPUS_PER_TASK $DEST/lineage.ms $DEST

