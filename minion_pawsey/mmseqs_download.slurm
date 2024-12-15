#!/bin/bash
#SBATCH --job-name=MMSeqsDld
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -o mmseqsdld-%j.out
#SBATCH -e mmseqsdld-%j.err

set -euo pipefail
eval "$(conda shell.bash hook)"
#conda activate atavide_lite
conda activate /scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite



DATE=`date +%Y%m%d`


TPD=$(mktemp -d -p .)
DB=UniRef50
mkdir -p $DB.$DATE/
mmseqs databases --threads 8 $DB $DB.$DATE/$DB $TPD
rm -rf $TPD

