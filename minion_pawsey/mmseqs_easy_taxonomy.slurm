#!/bin/bash

###############################################################
#                                                             #
# mmseqs easy taxonomy against UniRef50                       #
#                                                             #
# Note that you need a fasta file                             #
# Also make sure that the reads are labeled /1 /2             #
# (without a space)                                           #
#                                                             #
#                                                             #
###############################################################


#SBATCH --job-name=mmseqsLT
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH -o slurm_output/mmseqs_slurm/mmseqsET-%A_%a.out
#SBATCH -e slurm_output/mmseqs_slurm/mmseqsET-%A_%a.err


set -euo pipefail
eval "$(conda shell.bash hook)"
#conda activate atavide_lite
conda activate /scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite





# Output file is the base name. There will be several files with that name, so we write to $BGFS and then copy everything to the directory
# with that name once done

echo "Start " `date` >&2;

DB=UniRef50
DBSOURCE=UniRef50

if [[ ! -e $DBSOURCE ]]; then
	echo "Please edit $0 and correct the location of the $DB files" >&2;
	exit 1;
fi

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh

OUTDIR=mmseqs

R=$(head -n $SLURM_ARRAY_TASK_ID reads.txt | tail -n 1)
# convert fastq to fasta
# this throws an unbound error if FAFILEEND is not in DEFINITIONS
R1=${R/$FILEEND/$FAFILEEND}
OUTPUT=$OUTDIR/${R/$FILEEND/}
TMPOUTPUT=${R/$FILEEND/}

if [[ -e $OUTPUT ]]; then
	echo "$OUTPUT exists. Nothing to do" >&2;
	exit 0;
fi



mkdir -p $OUTPUT
echo "Running mmseqs easy-taxonomy fasta/$R1 $DBSOURCE/$DB $OUTPUT/$TMPOUTPUT tmp --threads 32" 2>&1;

TPD=$(mktemp -d -p .)
mmseqs easy-taxonomy fasta/$R1 $DBSOURCE/$DB $OUTPUT/$TMPOUTPUT $TPD --threads 32
rm -rf $TPD

find $OUTPUT/ -type f  | parallel -j 32 gzip

echo "Fin " `date` >&2;

