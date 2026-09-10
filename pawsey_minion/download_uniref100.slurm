#!/bin/bash
#SBATCH --job-name=MMSeqsDld100
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=7000M
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#SBATCH --partition=highmem

# NOTE: Setonix enforces MaxMemPerCPU (work=1840M, highmem=7900M). A bare --mem
# larger than cpus-per-task * MaxMemPerCPU does not fail - Slurm silently raises the
# CPU count to cover it, so the job queues for far more CPUs than it uses. Request
# memory with --mem-per-cpu so the two can never disagree.

# Use `mmseqs databases` to find information about the preformatted databases
# NOTE: Do not download the NT database unless you _really_ need it. The database is huge!


set -euo pipefail
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# test we have atavide_lite conda installed, and then activate it
eval "$(conda shell.bash hook)"
"$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh"
ATAVIDE_CONDA="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite"
conda activate $ATAVIDE_CONDA

# Derive thread count from the Slurm allocation rather than hardcoding it.
THREADS=${SLURM_CPUS_PER_TASK:-24}

CWD=`pwd`

#DB can also be  UniRef100 UniRef100 GTDB etc
DB=UniRef100


required_files=(
	UniRef100  UniRef100.dbtype  UniRef100_h  UniRef100_h.dbtype
	UniRef100_h.index  UniRef100.index  UniRef100.lookup  UniRef100_mapping
	UniRef100.source  UniRef100_taxonomy  UniRef100.version
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
mmseqs databases --threads "$THREADS" $DB $DESTINATION/$DB $TPD



