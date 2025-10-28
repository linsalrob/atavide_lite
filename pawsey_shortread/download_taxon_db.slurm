#!/bin/bash
#SBATCH --job-name=taxon_dld
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH -o slurm_output/taxon_download-%j.out
#SBATCH -e slurm_output/taxon_download-%j.err

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
