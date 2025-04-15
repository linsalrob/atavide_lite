#!/bin/bash
#SBATCH --job-name=taxon_dld
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH -o slurm_output/taxon_download-%j.out
#SBATCH -e slurm_output/taxon_download-%j.err

TAXDIR="/home/edwa0468/Databases/NCBI/taxonomy/current"
mkdir -p $TAXDIR
cd $TAXDIR
curl -LO https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar zxf taxdump.tar.gz && rm -f taxdump.tar.gz

