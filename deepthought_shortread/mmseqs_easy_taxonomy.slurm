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
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH -o slurm_output/mmseqs_slurm/mmseqsET-%A_%a.out
#SBATCH -e slurm_output/mmseqs_slurm/mmseqsET-%A_%a.err

eval "$(conda shell.bash hook)"
conda activate atavide_lite
set -euo pipefail

THREADS=${SLURM_CPUS_PER_TASK:-32}




# Output file is the base name. There will be several files with that name, so we write to $BGFS and then copy everything to the directory
# with that name once done

echo "Start " `date` >&2;

# Which MMseqs database to search. Override with e.g. MMSEQS_DB=UniRef100.
DB="${MMSEQS_DB:-UniRef50}"
DBSOURCE=/scratch/user/$USER/Databases/mmseqs/$DB

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
mkdir --parents $BGFS/output

R=$(head -n $SLURM_ARRAY_TASK_ID R1_reads.txt | tail -n 1)
# convert fastq to fasta
# this throws an unbound error if FAFILEEND is not in DEFINITIONS
R1=${R/$FILEEND/$FAFILEEND}
R2=${R1/_R1/_R2}
OUTPUT=$OUTDIR/${R/$FILEEND/}
TMPOUTPUT=${R/$FILEEND/}

# A sample is done only when its gzipped tophit report exists: the gzip runs after
# mmseqs finishes, so a cancelled or failed run never leaves one behind. Record which
# database produced the directory and refuse to reuse results from a different one.
if [[ -s "$OUTPUT/${TMPOUTPUT}_tophit_report.gz" ]]; then
	PREV_DB=$(cat "$OUTPUT/.mmseqs_db" 2>/dev/null || echo UniRef50)
	if [[ "$PREV_DB" == "$DB" ]]; then
		echo "$OUTPUT is complete (database $DB). Nothing to do" >&2;
		exit 0;
	fi
	echo "ERROR: $OUTPUT was built with database '$PREV_DB', not '$DB'." >&2;
	echo "       Remove it or choose a different output location before re-running." >&2;
	exit 3;
fi

mkdir -p $BGFS/$DB
rsync -a $DBSOURCE/ $BGFS/$DB/
cp fasta/$R1 fasta/$R2 $BGFS
echo "Copied " `date` >&2;

echo "Running mmseqs easy-taxonomy $BGFS/$R1 $BGFS/$R2 $BGFS/$DB/$DB $BGFS/output/$TMPOUTPUT tmp --threads $THREADS" >&2;

mmseqs easy-taxonomy $BGFS/$R1 $BGFS/$R2 $BGFS/$DB/$DB $BGFS/output/$TMPOUTPUT $(mktemp -d -p $BGFS) --threads $THREADS

find $BGFS/output -type f  | parallel -j $THREADS gzip

echo "Copying results from $BGFS/output* to $OUTPUT" >&2;

mkdir -p $OUTPUT
rsync -av $BGFS/output/ $OUTPUT/
echo "$DB" > "$OUTPUT/.mmseqs_db"

echo "Fin " `date` >&2;

