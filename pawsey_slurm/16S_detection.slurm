#!/bin/bash
#SBATCH --job-name=16S
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH -o slurm_output/sixteen_s/16S-%A_%a.out
#SBATCH -e slurm_output/sixteen_s/16S-%A_%a.err

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

# copy the 16S partie database

DBDIR=/scratch/$PAWSEY_ACCOUNT/$USER/Databases/partie
if [[ $SLURM_ARRAY_TASK_ID == 1 ]]; then
	if [ ! -e $DBDIR/copying_done ]; then
		mkdir -p $DBDIR
		rclone copy A18:dbs/partie/16SMicrobial.fasta.gz $DBDIR/
		touch $DBDIR/copying_done
	fi
else
	while [ ! -e $DBDIR/copying_done ]; do
		# sleep a random amount of time between 1 and 60 seconds. 
		# Hopefully this will also cause random offsets in reading the db
		sleep $(( (RANDOM % 30) + 1 ))
	done
fi

QC=fastq_fastp

R1=$(head -n $SLURM_ARRAY_TASK_ID R1_reads.txt | tail -n 1)
R2=${R1/_R1/_R2}

echo "Processing $R1 and $R2" >&2;

mkdir -p 16S/bam 16S/matches 16S/no_matches

O=${R1/$FILEEND/.bam}
if [[ $O == $R1 ]]; then
	echo "ERROR: OUTPUT $O and INPUT $R1 are the same. Are you sure the R1 files end $FILEEND?" >&2;
	exit 1;
fi

SIXS=$DBDIR/16SMicrobial.fasta.gz
if [[ ! -e $SIXS ]]; then
	echo "ERROR: $SIXS not found even though we have copied it?" >&2;
	exit 2
fi

echo -e "minimap2 -t 16 --split-prefix=tmp$$ -a -xsr $SIXS $QC/$R1 $QC/$R2 | samtools view -bh | samtools sort -o 16S/bam/$O -" >&2;
minimap2 -t 16 --split-prefix=tmp$$ -a -xsr $SIXS $QC/$R1 $QC/$R2 | samtools view -bh | samtools sort -o 16S/bam/$O -
samtools index 16S/bam/$O


# Reads that match the host
samtools fastq -F 3588 -f 65 16S/bam/$O | gzip -c > 16S/matches/$R1
samtools fastq -F 3588 -f 129 16S/bam/$O | gzip -c > 16S/matches/$R2


# Reads that do not match the host

samtools fastq -F 3584 -f 77 16S/bam/$O | gzip -c > 16S/no_matches/$R1
samtools fastq -F 3584 -f 141 16S/bam/$O | gzip -c > 16S/no_matches/$R2

# We don't currently save singleton reads, but these are they:
# samtools fastq -f 4 -F 1 $SAMPLE.bam > $SAMPLE.unmapped.single.fastq

samtools flagstats 16S/bam/$O

