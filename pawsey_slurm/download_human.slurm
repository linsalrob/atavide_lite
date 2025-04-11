#!/bin/bash
#SBATCH --job-name=human_download
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH -o slurm/human_download-%j.out
#SBATCH -e slurm/human_download-%j.err

DEST=/scratch/$PAWSEY_PROJECT/$USER/Databases/human

mkdir -p $DEST
cd $DEST


curl -LO https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
echo "Downloaded human db to $DEST"
