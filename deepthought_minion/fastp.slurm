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

# Same outputs as pawsey_minion/fastp.slurm: trimmed reads in fastq_fastp/ and the
# fastp JSON/HTML reports in fastp_output/.
mkdir -p fastq_fastp fastp_output
THREADS=${SLURM_CPUS_PER_TASK:-1}
echo "fastp -n 2 -l 50 -i fastq/$READS -o fastq_fastp/$READS -j fastp_output/$READS.json -h fastp_output/$READS.html --adapter_fasta ~/GitHubs/atavide_lite/adapters/ONT_RPB114.24.fna --thread $THREADS" >&2
fastp -n 2 -l 50 -i fastq/$READS -o fastq_fastp/$READS -j fastp_output/$READS.json -h fastp_output/$READS.html --adapter_fasta ~/GitHubs/atavide_lite/adapters/ONT_RPB114.24.fna --thread $THREADS
