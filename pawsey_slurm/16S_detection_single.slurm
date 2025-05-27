#!/bin/bash
#SBATCH --job-name=16S
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH -o slurm_output/sixteen_s/16S-single_%j.out
#SBATCH -e slurm_output/sixteen_s/16S-single_%j.err


# This script runs all the 16S in one go because its actually quite quick to do
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
SIXS=$DBDIR/16SMicrobial.fasta.gz
if [ ! -e $SIXS ]; then
	mkdir -p $DBDIR
	rclone copy A18:dbs/partie/16SMicrobial.fasta.gz $DBDIR/
fi

FQ=fastq
mkdir -p 16S/bam 

for R1 in $(cat R1_reads.txt); do
	R2=${R1/_R1/_R2}
	echo "Processing $R1 and $R2" >&2;


	O=${R1/$FILEEND/.bam}
	if [[ $O == $R1 ]]; then
		echo "ERROR: OUTPUT $O and INPUT $R1 are the same. Are you sure the R1 files end $FILEEND?" >&2;
		exit 1;
	fi

	echo -e "minimap2 -t 16 --split-prefix=tmp$$ -a -xsr $SIXS $FQ/$R1 $FQ/$R2 | samtools view -bh | samtools sort -o 16S/bam/$O -" >&2;
	minimap2 -t 16 --split-prefix=tmp$$ -a -xsr $SIXS $FQ/$R1 $FQ/$R2 | samtools view -bh | samtools sort -o 16S/bam/$O -
	samtools index 16S/bam/$O

	# make a json format output, one record per line
	echo -ne '{"'${R1/$FILEEND/}'": ' >> 16S/all_stats.json
	samtools flagstats -O json 16S/bam/$O | perl -pe 'chomp' >> 16S/all_stats.json
	echo "}" >> 16S/all_stats.json
done

python ~/GitHubs/atavide_lite/bin/summarise_16S.py --json 16S/all_stats.json --output 16S_counts.tsv
