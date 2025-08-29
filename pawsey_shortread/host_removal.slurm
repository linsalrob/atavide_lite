#!/bin/bash
#SBATCH --job-name=HostRemoval
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH -o slurm_output/host_slurm/host-%A_%a.out
#SBATCH -e slurm_output/host_slurm/host-%A_%a.err

set -euo pipefail

if [[ ! -n "$ATAVIDE_CONDA" ]]; then
    echo "Please define the ATAVIDE_CONDA environment variable with the location of your ATAVIDE_LITE conda installation" >&2;
    exit 2;
fi

eval "$(conda shell.bash hook)"
conda activate $ATAVIDE_CONDA

if [[ ! -e R1_reads.txt ]]; then
	echo "Please make a file with the R1 reads using this command:" >&2
	echo "find fastq -name \*R1\* > R1_reads.txt" >&2;
	exit 2;
fi

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED, HOSTFILE, etc" >&2
	exit 2;
fi

source DEFINITIONS.sh

QC=fastq_fastp

R1=$(head -n $SLURM_ARRAY_TASK_ID R1_reads.txt | tail -n 1)
R2=${R1/_R1/_R2}

echo "Processing $R1 and $R2" >&2;

mkdir -p $HOST $HOSTREMOVED host_bamfiles

O=${R1/$FILEEND/.bam}
if [[ $O == $R1 ]]; then
	echo "ERROR: OUTPUT $O and INPUT $R1 are the same. Are you sure the R1 files end $FILEEND?" >&2;
	exit 1;
fi

echo -e "minimap2 -t 16 --split-prefix=tmp$$ -a -xsr $HOSTFILE $QC/$R1 $QC/$R2 | samtools view -bh | samtools sort -o host_bamfiles/$O -" >&2;
minimap2 -t 16 --split-prefix=tmp$$ -a -xsr $HOSTFILE $QC/$R1 $QC/$R2 | samtools view -bh | samtools sort -o host_bamfiles/$O -
samtools index host_bamfiles/$O


# Reads that match the host
samtools fastq -F 3588 -f 65 host_bamfiles/$O | gzip -c > $HOST/$R1
samtools fastq -F 3588 -f 129 host_bamfiles/$O | gzip -c > $HOST/$R2


# Reads that do not match the host

samtools fastq -F 3584 -f 77 host_bamfiles/$O | gzip -c > $HOSTREMOVED/$R1
samtools fastq -F 3584 -f 141 host_bamfiles/$O | gzip -c > $HOSTREMOVED/$R2

# We don't currently save singleton reads, but these are they:
# samtools fastq -f 4 -F 1 $SAMPLE.bam > $SAMPLE.unmapped.single.fastq

