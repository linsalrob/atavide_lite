#!/usr/bin/env bash

# # Let's move the appropriate DEFINITION.sh
# Change this to the appropriate file

source DEFINITIONS.sh
# # Where are our metagenomic fastq.gz files?
DIR="$SOURCE"
# # Where is atavide_lite located on the HPC?
export PAWSEY_SRC=$HOME/GitHubs/atavide_lite/pawsey_shortread

#cp $PAWSEY_SRC/DEFINITIONS_human.sh DEFINITIONS.sh
if [[ ! -r "DEFINITIONS.sh" ]]; then
  printf 'ERROR: DEFINITIONS.sh is missing or unreadable in %s\n' "$PWD" >&2
  exit 1
fi

if [[ ! -d "$DIR" ]]; then
    printf 'ERROR: Input directory does not exist: %s\n' "$DIR" >&2
    exit 1
fi

if [[ ! -d "$PAWSEY_SRC" ]]; then
    printf 'ERROR: atavide_lite directory does not exist: %s\n' \
        "$PAWSEY_SRC" >&2
    exit 1
fi

: "${FILEEND:?FILEEND is not defined in DEFINITIONS.sh}"
: "${SOURCE:?SOURCE is not defined in DEFINITIONS.sh}"
: "${HOSTREMOVED:?HOSTREMOVED is not defined in DEFINITIONS.sh}"

printf 'Working directory: %s\n' "$PWD"
printf 'FASTQ directory:   %s\n' "$DIR"
printf 'SOURCE:            %s\n' "$SOURCE"
printf 'HOSTREMOVED:       %s\n' "$HOSTREMOVED"
printf 'FILEEND:           %s\n' "$FILEEND"

# # Make some files that will be useful for later
#find "$DIR" -maxdepth 1 -type f -name '*.fastq.gz' \
#  -printf '%f\n' > fastq_list.txt

#find "$DIR" -maxdepth 1 -type f -name '*.fastq.gz' \
#  -printf '%p\n' > fastq_list_full_path.txt

find "$DIR" -name \*_R1\* -printf "%f\n" > R1_reads.txt
export NUM_R1_READS=$(wc -l R1_reads.txt | cut -f 1 -d ' ')

if (( NUM_R1_READS == 0 )); then
    printf 'ERROR: No R1 files were found in %s\n' "$DIR" >&2
    exit 1
fi

echo "There are $NUM_R1_READS samples to process"

# Confirm that every listed R1 has a corresponding R2.
MISSING_PAIRS=0

while IFS= read -r R1; do
    R2=${R1/_R1/_R2}

    if [[ ! -s "$DIR/$R2" ]]; then
        printf 'ERROR: Missing or empty R2 for %s: %s\n' \
            "$R1" "$DIR/$R2" >&2
        MISSING_PAIRS=$((MISSING_PAIRS + 1))
    fi
done < R1_reads.txt

if (( MISSING_PAIRS > 0 )); then
    printf 'ERROR: Found %s incomplete read pairs\n' \
        "$MISSING_PAIRS" >&2
    exit 1
fi

printf 'All %s R1/R2 pairs were found\n' "$NUM_R1_READS"


# # Make some directories
RUN_DATE=$(date +%F)
FASTP_LOG_DIR="slurm_output/fastp_${RUN_DATE}_slurm"
MEGAHIT_LOG_DIR="slurm_output/megahit_${RUN_DATE}_slurm"
mkdir -p slurm_output/fastp_${RUN_DATE}_slurm 
mkdir -p slurm_output/megahit_${RUN_DATE}_slurm  



#VAMB_INSTALL=$(sbatch --parsable $PAWSEY_SRC/vamb_install.slurm)

# FASTP - Data quality control step
JOB=$(sbatch --parsable \
      --array="1-${NUM_R1_READS}:1" \
      --output="${FASTP_LOG_DIR}/fastp-%A_%a.out" \
      --error="${FASTP_LOG_DIR}/fastp-%A_%a.err" \
      "$PAWSEY_SRC/fastp.slurm"
)
printf 'FASTP array: %s\n' "$JOB"

# Host remove (if defined in the DEFINITION.sh file)
#HOSTJOB=$(sbatch --parsable --array=1-$NUM_R1_READS:1 --dependency=afterok:$JOB $PAWSEY_SRC/host_removal.slurm)

#MEGAHITHR=$(sbatch --parsable --dependency=afterok:$HOSTJOB $PAWSEY_SRC/megahit_allreads.slurm)

MEGAHITJOB=$(sbatch --parsable \
            --array="1-${NUM_R1_READS}:1" \
            --dependency="aftercorr:${JOB}" \
            --output="${MEGAHIT_LOG_DIR}/megahit-%A_%a.out" \
            --error="${MEGAHIT_LOG_DIR}/megahit-%A_%a.err" \
            "$PAWSEY_SRC/megahit.slurm"
)

printf 'MEGAHIT array: %s\n' "$MEGAHITJOB"