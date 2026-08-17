#!/usr/bin/env bash

# # Let's move the appropriate DEFINITION.sh
# Change this to the appropriate file

source DEFINITIONS.sh
# # Where are our metagenomic fastq.gz files?
DIR="$SOURCE"
# # Where is atavide_lite located on the HPC?
export PAWSEY_SRC=$HOME/GitHubs/atavide_lite/pawsey_shortread
export SOFTWARE_SRC=$MYPROJECT/software/miniconda3/
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
MMSEQS_LOG_DIR="slurm_output/mmseqs_${RUN_DATE}_slurm"
VAMB_LOG_DIR="slurm_output/vamb_${RUN_DATE}_slurm"
CHECKM_LOG="slurm_output/checkm_${RUN_DATE}"
 
#mkdir -p slurm_output/mmseqs_${RUN_DATE}_slurm
mkdir -p slurm_output/vamb_${RUN_DATE}_slurm 


#VAMB_INSTALL=$(sbatch --parsable \
#            --account="${PAWSEY_PROJECT}-gpu" \
#            $PAWSEY_SRC/vamb_install.slurm
#)


# Binning step #
#VCJOB=$(sbatch --parsable \
#            --dependency="afterok:${VAMB_INSTALL}" \
#            --output="${VAMB_LOG_DIR}/vamb_concat-%j.out" \
#            --error="${VAMB_LOG_DIR}/vamb_concat-%j.err" \
#            "$PAWSEY_SRC/vamb_concat.slurm"
#)

#VCJOB=$(sbatch --parsable \
#            --output="${VAMB_LOG_DIR}/vamb_concat-%j.out" \
#            --error="${VAMB_LOG_DIR}/vamb_concat-%j.err" \
#            "$PAWSEY_SRC/vamb_concat.slurm"
#)
printf 'VAM concatenation: %s\n' "$VCJOB"

#VMJOB=$(sbatch --parsable  \
#            --array="1-${NUM_R1_READS}:1" \
#            --dependency="afterok:${VCJOB}" \
#            --output="${VAMB_LOG_DIR}/map_reads-%A_%a.out" \
#            --error="${VAMB_LOG_DIR}/map_reads-%A_%a.err" \
#            "$PAWSEY_SRC/vamb_minimap.slurm"
#)

printf 'VAMB mapping array: %s\n' "$VMJOB"


#VAMBJOB=$(sbatch --parsable \
#            --dependency="afterok:${VMJOB}" \
#            --account="${PAWSEY_PROJECT}-gpu" \
#            --output="${VAMB_LOG_DIR}/vamb-%j.out" \
#            --error="${VAMB_LOG_DIR}/vamb-%j.err" \
#            "$PAWSEY_SRC/vamb.slurm"
#)


VAMBJOB=$(sbatch --parsable \
            --account="${PAWSEY_PROJECT}-gpu" \
            --output="${VAMB_LOG_DIR}/vamb-%j.out" \
            --error="${VAMB_LOG_DIR}/vamb-%j.err" \
            "$PAWSEY_SRC/vamb.slurm"
)

printf 'VAMB binning: %s\n' "$VAMBJOB"

CHECKMJOB=$(sbatch --parsable \
            --dependency="afterok:${VAMBJOB}" \
            --output="${CHECKM_LOG}_checkm-%j.out" \
            --error="${CHECKM_LOG}_checkm-%j.err" \
            "$PAWSEY_SRC/checkm.slurm" \
            "vamb/bins/" \
            "vamb/checkm"

)

printf 'CheckM: %s\n' "$CHECKMJOB"
printf '\nSubmission completed successfully.\n'
printf 'Dependency chain:\n'
