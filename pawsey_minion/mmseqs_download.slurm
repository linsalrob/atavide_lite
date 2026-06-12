#!/bin/bash
#SBATCH --job-name=MMSeqsDld
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -o mmseqsdld-%j.out
#SBATCH -e mmseqsdld-%j.err

set -euo pipefail
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# test we have atavide_lite conda installed, and then activate it
eval "$(conda shell.bash hook)"
"$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh"
ATAVIDE_CONDA="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite"
conda activate $ATAVIDE_CONDA

DATE=`date +%Y%m%d`

TPD=$(mktemp -d -p /scratch/$PAWSEY_PROJECT/$USER/tmp)
DB=UniRef50
mkdir -p $DB.$DATE/
mmseqs databases --threads 8 $DB $DB.$DATE/$DB $TPD
rm -rf $TPD

