#!/bin/bash
#SBATCH --job-name=MMSeqsDld
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -o mmseqsdld-%j.out
#SBATCH -e mmseqsdld-%j.err

set -euo pipefail
if [[ ! -n "$ATAVIDE_CONDA" ]]; then
    echo "Please include the ATAVIDE_CONDA location on your slurm submission: --export=ATAVIDE_CONDA=$ATAVIDE_CONDA" >&2
    exit 2;
fi
eval "$(conda shell.bash hook)"
conda activate $ATAVIDE_CONDA

DATE=`date +%Y%m%d`


TPD=$(mktemp -d -p /scratch/$PAWSEY_PROJECT/$USER/tmp)
DB=UniRef50
mkdir -p $DB.$DATE/
mmseqs databases --threads 8 $DB $DB.$DATE/$DB $TPD
rm -rf $TPD

