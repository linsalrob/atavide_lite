#!/bin/bash

###############################################################
#                                                             #
# mmseqs easy taxonomy against UniRef100 (MMSEQS_DB)         #
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


set -euo pipefail
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# test we have atavide_lite conda installed, and then activate it
eval "$(conda shell.bash hook)"
"$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh"
ATAVIDE_CONDA="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite"
conda activate $ATAVIDE_CONDA

echo "Start " `date` >&2;

# Which MMseqs database to search. Override with e.g. MMSEQS_DB=UniRef50 (smaller and faster).
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

R=$(head -n $SLURM_ARRAY_TASK_ID R1_reads.txt | tail -n 1)

TMPDIR=/scratch/$PAWSEY_PROJECT/$USER/tmp
OUTDIR=mmseqs/${R/$FILEEND/}
mkdir --parents $OUTDIR $TMPDIR

# convert fastq to fasta
FAFILEEND=${FILEEND/fastq/fasta}
R1=${R/$FILEEND/$FAFILEEND}
R2=${R1/_R1/_R2}
OUTPUT=$OUTDIR/${R/$FILEEND/}
TMPOUTPUT=${R/$FILEEND/}

# A sample is done only when its gzipped tophit report exists: the gzip runs after
# mmseqs finishes, so a cancelled or failed run never leaves one behind. Record which
# database produced the directory and refuse to reuse results from a different one.
if [[ -s "${OUTPUT}_tophit_report.gz" ]]; then
	# directories from before the label existed were built with the old default, UniRef50
	PREV_DB=$(cat "$OUTDIR/.mmseqs_db" 2>/dev/null || echo UniRef50)
	if [[ "$PREV_DB" == "$DB" ]]; then
		echo "$OUTPUT is complete (database $DB). Nothing to do" >&2;
		exit 0;
	fi
	echo "ERROR: $OUTDIR was built with database '$PREV_DB', not '$DB'." >&2;
	echo "       To keep using it, resubmit with --export=ALL,MMSEQS_DB=$PREV_DB; otherwise" >&2;
	echo "       remove it or choose a different output location before re-running." >&2;
	exit 3;
fi
echo "$DB" > "$OUTDIR/.mmseqs_db"

# cp fasta/$R1 fasta/$R2 $BGFS

# MMseqs defaults --split-memory-limit to all the memory in /proc/meminfo, i.e. the whole
# node rather than this job's allocation, and can then be OOM-killed in prefilter against
# UniRef100. Cap it at 80% of the allocation, as pawsey_minion does.
if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
	SPLIT_MEM_ARG=(--split-memory-limit "$(( SLURM_MEM_PER_NODE * 4 / 5 ))M")
elif [[ -n "${SLURM_MEM_PER_CPU:-}" && -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
	SPLIT_MEM_ARG=(--split-memory-limit "$(( SLURM_MEM_PER_CPU * SLURM_CPUS_PER_TASK * 4 / 5 ))M")
else
	SPLIT_MEM_ARG=()
	echo "WARNING: could not determine allocated memory; mmseqs will assume all node memory" >&2
fi

# create the temporary directory once, so the path we log is the one mmseqs uses
TPD=$(mktemp -d -p $TMPDIR)
echo "Running  mmseqs easy-taxonomy fasta/$R1 fasta/$R2 $DBSOURCE/$DB $OUTPUT $TPD --threads 64 ${SPLIT_MEM_ARG[*]}" >&2;
mmseqs easy-taxonomy fasta/$R1 fasta/$R2 $DBSOURCE/$DB $OUTPUT $TPD --threads 64 "${SPLIT_MEM_ARG[@]}"

# leave the database label readable, and overwrite stale .gz files left by a cancelled attempt
find $OUTDIR -type f ! -name '*.gz' ! -name .mmseqs_db | parallel -j 32 gzip -f

echo "Fin " `date` >&2;

