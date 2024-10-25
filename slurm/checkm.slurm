#!/bin/bash
#SBATCH --job-name=checkm
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH -o slurm_output/checkm-%j.out
#SBATCH -e slurm_output/checkm-%j.err

if [[ -z $CHECKM_DATA_PATH ]]; then 
	echo "Please define CHECKM_DATA_PATH with your databases" >&2; 
	exit 1;
fi

if [[ $# -lt 2 ]]; then
	echo "$0 <source directory> <destination directory>" >&2
	exit 1
fi

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate checkm

# this %/ removes the trailing / if there is one
SOURCE=${1%/}
DEST=${2%/}
LOCALSOURCE=`basename $SOURCE`

echo "Copying source to $BGFS" >&2
rsync -az $SOURCE/ $BGFS/$LOCALSOURCE/

echo "Running checkm on $BGFS/$LOCALSOURCE and output in $BGFS/checkm"  >&2;
checkm lineage_wf -t 16 $BGFS/$LOCALSOURCE $BGFS/checkm
checkm qa -t 16 $BGFS/checkm/lineage.ms $BGFS/checkm
rsync -az $BGFS/checkm/ $DEST/

