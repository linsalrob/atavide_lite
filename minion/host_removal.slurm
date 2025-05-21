#!/bin/bash
#SBATCH --job-name=HostRemoval
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH -o slurm_output/host_slurm/host-%A_%a.out
#SBATCH -e slurm_output/host_slurm/host-%A_%a.err

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate atavide_lite


if [[ ! -e reads.txt ]]; then
	echo "Please make a file with the reads using this command:" >&2
	echo "find fastq -type f -printf "%f\n" > reads.txt" >&2;
	exit 2;
fi

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED, HOSTFILE, etc" >&2
	exit 2;
fi

source DEFINITIONS.sh


QC=fastq_fastp
READS=$(head -n $SLURM_ARRAY_TASK_ID reads.txt | tail -n 1)

echo "Processing $READS" >&2;

mkdir -p $HOST $HOSTREMOVED

O=${READS/$FILEEND/.bam}
if [[ $O == $READS ]]; then
	echo "ERROR: OUTPUT $O and INPUT $READS are the same. Are you sure the files end $FILEEND?" >&2;
	exit 1;
fi
TMPDIR=$(mktemp -d -p /scratch/$PAWSEY_PROJECT/$USER/tmp)
echo -e "minimap2 -t 16 --split-prefix=$TMPDIR/tmp$$ -a -x map-ont $HOSTFILE fastq_fastp/$READS | samtools view -bh | samtools sort -o $HOST/$O -\n" >&2;
minimap2 -t 64 --secondary=no --split-prefix=$TMPDIR/tmp$$ -a -x map-ont $HOSTFILE fastq_fastp/$READS | samtools view -bh | samtools sort -o $HOST/$O -
samtools index $HOST/$O


# Reads that match the host
samtools fastq -F 3588 $HOST/$O | gzip -c > $HOST/$READS


# Reads that do not match the host

samtools fastq -F 3584 -f 4 $HOST/$O | gzip -c > $HOSTREMOVED/$READS

