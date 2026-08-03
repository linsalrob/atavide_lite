#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH -o slurm_output/fastp_slurm/fastp-%A_%a.out
#SBATCH -e slurm_output/fastp_slurm/fastp-%A_%a.err

set -euo pipefail

eval "$(conda shell.bash hook)"
conda activate atavide_lite


if [[ ! -e reads.txt ]]; then 
	echo "Please make a file with the reads using this command:" >&2
	echo "find fastq -type f -printf "%f\n" > reads.txt" >&2;
	exit 2;
fi

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh

READS=$(head -n $SLURM_ARRAY_TASK_ID reads.txt | tail -n 1)

mkdir -p fastq_fastp fastq_fastp_output
echo "fastp -n 1 -l 100 -i fastq/$READS -o output/$READS --adapter_fasta ~/atavide_lite/adapters/ONT_Rapid_Barcodes.fna " >&2
fastp -n 2 -l 50 -i fastq/$READS -o fastq_fastp/$READS -h fastq_fastp_output/${READS/$FILEEND/}.html -j fastq_fastp_output/${READS/$FILEEND/}.json --adapter_fasta ~/atavide_lite/adapters/ONT_RPB114.24.fna

pigz fastq_fastp_output/${READS/$FILEEND/}.html &
pigz fastq_fastp_output/${READS/$FILEEND/}.json &

wait
