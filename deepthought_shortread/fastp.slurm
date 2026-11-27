#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -o slurm_output/fastp_slurm/fastp-%A_%a.out
#SBATCH -e slurm_output/fastp_slurm/fastp-%A_%a.err

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate atavide_lite


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

# Same outputs as pawsey_shortread/fastp.slurm: trimmed reads in fastq_fastp/ and
# the fastp JSON/HTML reports in fastp_output/.
mkdir -p fastq_fastp fastp_output

if [[ -s "fastq_fastp/$R1" && \
      -s "fastq_fastp/$R2" && \
      -s "fastp_output/$R1.json" && \
      -s "fastp_output/$R1.html" ]]; then
    echo "fastp outputs already exist for $R1. Skipping." >&2
    exit 0
fi

THREADS=${SLURM_CPUS_PER_TASK:-4}

mkdir -p $BGFS/output
mkdir -p $BGFS/fastq
cp ~/GitHubs/atavide_lite/adapters/IlluminaAdapters.fa $BGFS

cp $SOURCE/$R1 $SOURCE/$R2 $BGFS/fastq

echo "fastp -n 1 -l 100 -i $BGFS/fastq/$R1 -I $BGFS/fastq/$R2 -o $BGFS/output/$R1 -O $BGFS/output/$R2 -j $BGFS/output/$R1.json -h $BGFS/output/$R1.html --adapter_fasta $BGFS/IlluminaAdapters.fa --thread $THREADS" >&2
fastp -n 1 -l 100 -i $BGFS/fastq/$R1 -I $BGFS/fastq/$R2 -o $BGFS/output/$R1 -O $BGFS/output/$R2 -j $BGFS/output/$R1.json -h $BGFS/output/$R1.html --adapter_fasta $BGFS/IlluminaAdapters.fa --thread $THREADS

mv $BGFS/output/$R1.json $BGFS/output/$R1.html fastp_output
mv $BGFS/output/$R1 $BGFS/output/$R2  fastq_fastp
