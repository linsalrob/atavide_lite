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
#SBATCH --cpus-per-task=1
#SBATCH --mem=220G
#SBATCH -o slurm_output/mmseqs_slurm/mmseqsET-%A_%a.out
#SBATCH -e slurm_output/mmseqs_slurm/mmseqsET-%A_%a.err

eval "$(conda shell.bash hook)"

if [[ ! -n "$ATAVIDE_CONDA" ]]; then
    echo "Please define the ATAVIDE_CONDA environment variable with the location of your ATAVIDE_LITE conda installation" >&2;
    exit 2;
fi

eval "$(conda shell.bash hook)"
conda activate $ATAVIDE_CONDA

echo "Start " `date` >&2;

DB=UniRef50
DBSOURCE=/scratch/$PAWSEY_PROJECT/$USER/Databases/$DB
TMPDIR=/scratch/$PAWSEY_PROJECT/$USER/tmp

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
mkdir --parents $OUTDIR $TMPDIR

R=$(head -n $SLURM_ARRAY_TASK_ID R1_reads.txt | tail -n 1)
# convert fastq to fasta
# this throws an unbound error if FAFILEEND is not in DEFINITIONS
R1=${R/$FILEEND/$FAFILEEND}
R2=${R1/_R1/_R2}
OUTPUT=$OUTDIR/${R/$FILEEND/}
TMPOUTPUT=${R/$FILEEND/}

if [[ -e $OUTPUT ]]; then
	echo "$OUTPUT exists. Nothing to do" >&2;
	exit 0;
fi

# cp fasta/$R1 fasta/$R2 $BGFS

echo "Running  mmseqs easy-taxonomy fasta/$R1 fasta/$R2 $DBSOURCE/$DB $OUTPUT $(mktemp -d -p $TMPDIR) --threads 64" >&2;
mmseqs easy-taxonomy fasta/$R1 fasta/$R2 $DBSOURCE/$DB $OUTPUT $(mktemp -d -p $TMPDIR) --threads 64

find $OUTPUT -type f  | parallel -j 32 gzip

echo "Fin " `date` >&2;

