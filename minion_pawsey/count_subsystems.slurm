#!/bin/bash
#SBATCH --job-name=CountSS
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH -o slurm_output/count_ss-%j.out
#SBATCH -e slurm_output/count_ss-%j.err

perl ~/atavide_lite/bin/count_subsystems.pl -d mmseqs
