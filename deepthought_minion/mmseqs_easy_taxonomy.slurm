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
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH -o slurm_output/mmseqs_slurm/mmseqsET-%A_%a.out
#SBATCH -e slurm_output/mmseqs_slurm/mmseqsET-%A_%a.err


set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate atavide_lite

THREADS=${SLURM_CPUS_PER_TASK:-32}

echo "Start " `date` >&2;

# Which MMseqs database to search. Override with e.g. MMSEQS_DB=UniRef50 (smaller and faster).
DB="${MMSEQS_DB:-UniRef100}"
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
FAFILEEND=${FILEEND/.fastq/.fasta}
R1=${R/$FILEEND/$FAFILEEND}
OUTPUT=$OUTDIR/${R/$FILEEND/}
TMPOUTPUT=${R/$FILEEND/}

# A sample is done only when its gzipped tophit report exists: the gzip runs after
# mmseqs finishes, so a cancelled or failed run never leaves one behind. Record which
# database produced the directory and refuse to reuse results from a different one.
if [[ -s "$OUTPUT/${TMPOUTPUT}_tophit_report.gz" ]]; then
	# directories from before the label existed were built with the old default, UniRef50
	PREV_DB=$(cat "$OUTPUT/.mmseqs_db" 2>/dev/null || echo UniRef50)
	if [[ "$PREV_DB" == "$DB" ]]; then
		echo "$OUTPUT is complete (database $DB). Nothing to do" >&2;
		exit 0;
	fi
	echo "ERROR: $OUTPUT was built with database '$PREV_DB', not '$DB'." >&2;
	echo "       To keep using it, resubmit with --export=ALL,MMSEQS_DB=$PREV_DB; otherwise" >&2;
	echo "       remove it or choose a different output location before re-running." >&2;
	exit 3;
fi


mkdir -p $BGFS/$OUTPUT 

TPD=$(mktemp -d -p $BGFS)
echo "Running mmseqs easy-taxonomy fasta/$R1 $BGFS/$DB/$DB $BGFS/$OUTPUT/$TMPOUTPUT $TPD --threads $THREADS" >&2;
mmseqs easy-taxonomy fasta/$R1 $BGFS/$DB/$DB $BGFS/$OUTPUT/$TMPOUTPUT $TPD --threads $THREADS

find $BGFS/$OUTPUT/ -type f  | parallel -j $THREADS gzip
rsync -a $BGFS/$OUTPUT/ $OUTPUT/
echo "$DB" > "$OUTPUT/.mmseqs_db"

echo "Fin " `date` >&2;

