#!/bin/bash

###############################################################
#                                                             #
# mmseqs easy taxonomy against UniRef100                       #
#                                                             #
# Note that you need a fasta file                             #
# Also make sure that the reads are labeled /1 /2             #
# (without a space)                                           #
#                                                             #
#                                                             #
###############################################################


#SBATCH --job-name=mmseqsET100
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=480G
#SBATCH --partition=highmem
#SBATCH -o slurm_output/mmseqs_slurm/%x-%A_%a.out
#SBATCH -e slurm_output/mmseqs_slurm/%x-%A_%a.err


set -euo pipefail
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# test we have atavide_lite conda installed, and then activate it
eval "$(conda shell.bash hook)"
"$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh"
ATAVIDE_CONDA="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite"
conda activate $ATAVIDE_CONDA

echo "Start " `date` >&2;

DB=UniRef100
DBSOURCE=/scratch/$PAWSEY_PROJECT/$USER/Databases/$DB

if [[ ! -e $DBSOURCE ]]; then
	echo "Please edit $0 and correct the location of the $DB files. We expected them in $DBSOURCE" >&2;
	exit 1;
fi

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh

R=$(head -n $SLURM_ARRAY_TASK_ID R1_reads.txt | tail -n 1)

TMPDIR=/scratch/$PAWSEY_PROJECT/$USER/tmp
OUTDIR=mmseqs_uniref100/mmseqs/${R/$FILEEND/}
mkdir --parents $OUTDIR $TMPDIR

# convert fastq to fasta
FAFILEEND=${FILEEND/fastq/fasta}
R1=${R/$FILEEND/$FAFILEEND}
R2=${R1/_R1/_R2}
OUTPUT=$OUTDIR/${R/$FILEEND/}
TMPOUTPUT=${R/$FILEEND/}

if [[ -e $OUTPUT ]]; then
	echo "$OUTPUT exists. Nothing to do" >&2;
	exit 0;
fi

# cp fasta/$R1 fasta/$R2 $BGFS

echo "Running  mmseqs easy-taxonomy fasta/$R1 fasta/$R2 $DBSOURCE/$DB $OUTPUT $(mktemp -d -p $TMPDIR) --threads $SLURM_CPUS_PER_TASK" >&2;
mmseqs easy-taxonomy fasta/$R1 fasta/$R2 $DBSOURCE/$DB $OUTPUT $(mktemp -d -p $TMPDIR) --threads $SLURM_CPUS_PER_TASK

find $OUTDIR -type f  | parallel -j $SLURM_CPUS_PER_TASK gzip

echo "Fin " `date` >&2;

