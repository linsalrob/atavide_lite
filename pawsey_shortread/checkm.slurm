#!/bin/bash
#SBATCH --job-name=checkm
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB
#SBATCH --partition=highmem



set -euo pipefail
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# test we have atavide_lite conda installed, and then activate it
eval "$(conda shell.bash hook)"
"$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh"
ATAVIDE_CONDA="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite"
conda activate $ATAVIDE_CONDA

if [[ -z ${CHECKM_DATA_PATH:-} ]]; then 
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

# this %/ removes the trailing / if there is one
SOURCE=${1%/}
DEST=${2%/}
LOCALSOURCE=`basename $SOURCE`


echo "Running checkm on $SOURCE and output in checkm"  >&2;
checkm lineage_wf -t $SLURM_CPUS_PER_TASK $SOURCE $DEST
checkm qa -t $SLURM_CPUS_PER_TASK $DEST/lineage.ms $DEST

