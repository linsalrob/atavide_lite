#!/bin/bash
#SBATCH --job-name=MMSeqsDld
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH -o mmseqsdld-%j.out
#SBATCH -e mmseqsdld-%j.err

# Use `mmseqs databases` to find information about the preformatted databases
# NOTE: Do not download the NT database unless you _really_ need it. The database is huge!

set -euo pipefail

if [[ ! -n "$ATAVIDE_CONDA" ]]; then
    echo "Please define the ATAVIDE_CONDA environment variable with the location of your ATAVIDE_LITE conda installation" >&2;
    exit 2;
fi

eval "$(conda shell.bash hook)"
conda activate $ATAVIDE_CONDA


#DB can also be  UniRef50 UniRef100 GTDB etc
DB=UniRef50

DESTINATION=/scratch/$PAWSEY_PROJECT/$USER/Databases/$DB
TMP=/scratch/$PAWSEY_PROJECT/$USER/tmp
mkdir -p $DESTINATION $TMP

TPD=$(mktemp -d -p $TMP)
echo -e "Downloading to $DESTINATION/$DB and using temp location $TPD.\nNote we download the tarball and then uncompress, so there is usually nothing in the output for a while!" >&2
mmseqs databases --threads 8 $DB $DESTINATION/$DB $TPD



