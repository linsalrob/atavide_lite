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
CWD=`pwd`

#DB can also be  UniRef50 UniRef100 GTDB etc
DB=UniRef50


required_files=(
	UniRef50  UniRef50.dbtype  UniRef50_h  UniRef50_h.dbtype
	UniRef50_h.index  UniRef50.index  UniRef50.lookup  UniRef50_mapping
	UniRef50.source  UniRef50_taxonomy  UniRef50.version
)

missing=0

DESTINATION=/scratch/$PAWSEY_PROJECT/$USER/Databases/$DB
TMP=/scratch/$PAWSEY_PROJECT/$USER/tmp
mkdir -p $DESTINATION $TMP

cd $DESTINATION

for file in "${required_files[@]}"; do
  if [[ ! -f "$DB/$file" ]]; then
    echo "Missing: $DB/$file"
    missing=1
  fi
done

if [[ $missing -eq 0 ]]; then
  echo "All required files are present."
  exit 0
fi

echo "Downloading required files"

cd $CWD
# delete the old Destination and make the dir to remove partial files.

rm -rf $DESTINATION && mkdir -p $DESTINATION
TPD=$(mktemp -d -p $TMP)
echo -e "Downloading to $DESTINATION/$DB and using temp location $TPD.\nNote we download the tarball and then uncompress, so there is usually nothing in the output for a while!" >&2
mmseqs databases --threads 8 $DB $DESTINATION/$DB $TPD



