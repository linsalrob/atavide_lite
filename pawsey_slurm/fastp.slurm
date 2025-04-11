#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -o slurm_output/fastp_slurm/fastp-%A_%a.out
#SBATCH -e slurm_output/fastp_slurm/fastp-%A_%a.err


set -euo pipefail


if [[ ! -n "$ATAVIDE_CONDA" ]]; then
    echo "Please define the ATAVIDE_CONDA environment variable with the location of your ATAVIDE_LITE conda installation" >&2;
    exit 2;
fi



eval "$(conda shell.bash hook)"
conda activate $ATAVIDE_CONDA


if [[ ! -e R1_reads.txt ]]; then 
	echo "Please make a file with the R1 reads using this command:" >&2
	echo "find fastq -name \*R1\* -printf "%f\n" > R1_reads.txt" >&2;
	exit 2;
fi

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh

cp ~/atavide_lite/adapters/IlluminaAdapters.fa $BGFS

R1=$(head -n $SLURM_ARRAY_TASK_ID R1_reads.txt | tail -n 1)
R2=${R1/_R1/_R2}

echo "fastp -n 1 -l 100 -i fastq/$R1 -I fastq/$R2 -o fastq_fastp/$R1 -O fastq_fastp/$R2 --adapter_fasta $HOME/atavide_lite/adapters//IlluminaAdapters.fa " >&2
fastp -n 1 -l 100 -i fastq/$R1 -I fastq/$R2 -o fastq_fastp/$R1 -O fastq_fastp/$R2 --adapter_fasta $HOME/atavide_lite/adapters/IlluminaAdapters.fa 

