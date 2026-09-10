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
#SBATCH --mem-per-cpu=1800M
#SBATCH -o slurm_output/mmseqs_slurm/mmseqsET-%A_%a.out
#SBATCH -e slurm_output/mmseqs_slurm/mmseqsET-%A_%a.err

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

# MMseqs defaults --split-memory-limit to "all available system memory", which it reads
# from /proc/meminfo - i.e. the whole node (256 GB on a Setonix work node), not the memory
# actually granted to this job. It therefore decides the target DB needs no splitting and
# gets OOM-killed during prefilter. Derive the real limit from the allocation and leave
# some headroom for the rest of the pipeline.
if [[ -n "${SLURM_MEM_PER_CPU:-}" ]]; then
	ALLOC_MB=$(( SLURM_MEM_PER_CPU * THREADS ))
elif [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
	ALLOC_MB=$SLURM_MEM_PER_NODE
else
	ALLOC_MB=0
fi
if [[ "$ALLOC_MB" -gt 0 ]]; then
	SPLIT_MEM_MB=$(( ALLOC_MB * 4 / 5 ))
	SPLIT_MEM_ARG=(--split-memory-limit "${SPLIT_MEM_MB}M")
	echo "Allocation ${ALLOC_MB}M; capping mmseqs --split-memory-limit at ${SPLIT_MEM_MB}M" >&2
else
	SPLIT_MEM_ARG=()
	echo "WARNING: could not determine allocated memory; mmseqs will assume all node memory" >&2
fi

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
echo "Running mmseqs easy-taxonomy fasta/$R1 $DBSOURCE/$DB $OUTPUT/$TMPOUTPUT $TPD --threads $THREADS ${SPLIT_MEM_ARG[*]}" 2>&1;
mmseqs easy-taxonomy fasta/$R1 $DBSOURCE/$DB $OUTPUT/$TMPOUTPUT $TPD --threads "$THREADS" "${SPLIT_MEM_ARG[@]}"
rm -rf $TPD

find $OUTPUT/ -type f  | parallel -j 32 gzip

echo "Fin " `date` >&2;

