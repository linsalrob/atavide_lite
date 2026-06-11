#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH -o slurm_output/fastp_slurm/fastp-%A_%a.out
#SBATCH -e slurm_output/fastp_slurm/fastp-%A_%a.err

set -euo pipefail


set -euo pipefail
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# test we have atavide_lite conda installed, and then activate it
eval "$(conda shell.bash hook)"
"$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh"
ATAVIDE_CONDA="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite"
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


R1=$(head -n $SLURM_ARRAY_TASK_ID R1_reads.txt | tail -n 1)
R2=${R1/_R1/_R2}

mkdir -p fastq_fastp fastp_output
echo "fastp -n 1 -l 100 -i fastq/$R1 -I fastq/$R2 -o fastq_fastp/$R1 -O fastq_fastp/$R2 --adapter_fasta $HOME/atavide_lite/adapters//IlluminaAdapters.fa --thread 16" >&2
fastp -n 1 -l 100 -i fastq/$R1 -I fastq/$R2 -o fastq_fastp/$R1 -O fastq_fastp/$R2 -j fastp_output/$R1.json -h fastp_output/$R1.html --adapter_fasta $HOME/atavide_lite/adapters/IlluminaAdapters.fa --thread 16

