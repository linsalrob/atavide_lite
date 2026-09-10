#!/bin/bash

###############################################################
#                                                             #
# mmseqs easy taxonomy against UniRef100                      #
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
#SBATCH --mem-per-cpu=7000M
#SBATCH -o slurm_output/mmseqs_slurm/mmseqsET-%A_%a.out
#SBATCH -e slurm_output/mmseqs_slurm/mmseqsET-%A_%a.err
#SBATCH --partition=highmem

# NOTE: Setonix enforces MaxMemPerCPU (work=1840M, highmem=7900M). A bare --mem
# larger than cpus-per-task * MaxMemPerCPU does not fail - Slurm silently raises the
# CPU count to cover it, so the job queues for far more CPUs than it uses. Request
# memory with --mem-per-cpu so the two can never disagree.

set -euo pipefail
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# test we have atavide_lite conda installed, and then activate it
eval "$(conda shell.bash hook)"
"$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh"
ATAVIDE_CONDA="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite"
conda activate $ATAVIDE_CONDA

# Derive thread count from the Slurm allocation rather than hardcoding it.
THREADS=${SLURM_CPUS_PER_TASK:-32}

echo "Start " `date` >&2;

# Which MMseqs database to search. Override with e.g. MMSEQS_DB=UniRef50 (smaller and
# faster, so more samples run in parallel). Default is unchanged.
DB="${MMSEQS_DB:-UniRef100}"
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

OUTDIR=mmseqs

R=$(head -n $SLURM_ARRAY_TASK_ID reads.txt | tail -n 1)
# convert fastq to fasta
FAFILEEND=${FILEEND/.fastq/.fasta}
R1=${R/$FILEEND/$FAFILEEND}
OUTPUT=$OUTDIR/${R/$FILEEND/}
TMPOUTPUT=${R/$FILEEND/}

if [[ -e $OUTPUT ]]; then
	echo "$OUTPUT exists. Nothing to do" >&2;
	exit 0;
fi



mkdir -p $OUTPUT

TPD=$(mktemp -d -p /scratch/$PAWSEY_PROJECT/$USER/tmp)
echo "Running mmseqs easy-taxonomy fasta/$R1 $DBSOURCE/$DB $OUTPUT/$TMPOUTPUT $TPD --threads $THREADS" 2>&1;
mmseqs easy-taxonomy fasta/$R1 $DBSOURCE/$DB $OUTPUT/$TMPOUTPUT $TPD --threads "$THREADS"
rm -rf $TPD

find $OUTPUT/ -type f  | parallel -j 32 gzip

echo "Fin " `date` >&2;

