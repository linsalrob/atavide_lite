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
conda activate atavide_lite

echo "Start " `date` >&2;

DB=UniRef50
DBSOURCE=$HOME/Databases/mmseqs/$DB

if [[ ! -e $DBSOURCE ]]; then
	echo "Please edit $0 and correct the location of the $DB files. We expected them in $DBSOURCE" >&2;
	exit 1;
fi

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh

rsync -a $DBSOURCE/ $BGFS/$DB/
OUTDIR=mmseqs
mkdir -p $OUTDIR

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


mkdir -p $BGFS/$OUTPUT 

TPD=$(mktemp -d -p $BGFS)
echo "Running mmseqs easy-taxonomy fasta/$R1 $BGFS/$DB/$DB $BGFS/$OUTPUT/$TMPOUTPUT $TPD --threads 32" 2>&1;
mmseqs easy-taxonomy fasta/$R1 $BGFS/$DB/$DB $BGFS/$OUTPUT/$TMPOUTPUT $TPD --threads 32

find $BGFS/$OUTPUT/ -type f  | parallel -j 32 gzip
rsync -a $BGFS/$OUTPUT/ $OUTPUT/

echo "Fin " `date` >&2;

