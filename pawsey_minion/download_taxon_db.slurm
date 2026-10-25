#!/bin/bash
#SBATCH --job-name=taxon_dld
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1800M
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err

# NOTE: Setonix enforces MaxMemPerCPU (work=1840M, highmem=7900M). A bare --mem
# larger than cpus-per-task * MaxMemPerCPU does not fail - Slurm silently raises the
# CPU count to cover it, so the job queues for far more CPUs than it uses. Request
# memory with --mem-per-cpu so the two can never disagree.

TAXDIR="/scratch/$PAWSEY_PROJECT/$USER/Databases/NCBI/taxonomy/current"

mkdir -p $TAXDIR
cd $TAXDIR;

required_files=(
  citations.dmp  delnodes.dmp  division.dmp  gc.prt
  gencode.dmp  images.dmp  merged.dmp  names.dmp
  nodes.dmp  readme.txt
)

missing=0

for file in "${required_files[@]}"; do
  if [[ ! -f "$file" ]]; then
    echo "Missing: $file"
    missing=1
  fi
done

if [[ $missing -eq 1 ]]; then
  echo "One or more files are missing. Downloading taxdump.tar.gz..."
  curl -LO https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
  tar zxf taxdump.tar.gz && rm -f taxdump.tar.gz
else
  echo "All required files are present."
fi

echo "Done"
