#!/bin/bash
#SBATCH --job-name=checkm
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --mem-per-cpu=1800M
#SBATCH -o slurm_output/checkm-%j.out
#SBATCH -e slurm_output/checkm-%j.err

# NOTE: Setonix enforces MaxMemPerCPU (work=1840M). A bare --mem larger than
# cpus-per-task * MaxMemPerCPU does not fail - Slurm silently raises the CPU count to
# cover it, which also changes $SLURM_CPUS_PER_TASK and therefore THREADS below.


set -euo pipefail
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# test we have atavide_lite conda installed, and then activate it
eval "$(conda shell.bash hook)"
"$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh"
ATAVIDE_CONDA="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite"
conda activate $ATAVIDE_CONDA

if [[ -z "${CHECKM_DATA_PATH:-}" ]]; then 
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
THREADS=${SLURM_CPUS_PER_TASK:-16}


echo "Running checkm on $LOCALSOURCE and output in $DEST"  >&2;
checkm lineage_wf -t "$THREADS" "$SOURCE" "$DEST/checkm"
checkm qa -t "$THREADS" "$DEST/checkm/lineage.ms" "$DEST/checkm"

