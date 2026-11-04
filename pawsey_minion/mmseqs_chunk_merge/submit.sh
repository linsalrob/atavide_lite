#!/bin/bash
# Submit the chunked MMseqs workflow from an atavide analysis directory.

set -euo pipefail
[[ -s DEFINITIONS.sh && -s reads.txt ]] || {
    echo 'ERROR: run from an analysis directory containing DEFINITIONS.sh and reads.txt.' >&2
    exit 2
}
# shellcheck source=/dev/null
source DEFINITIONS.sh
: "${FILEEND:?DEFINITIONS.sh must set FILEEND.}"

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORK_ROOT=${WORK_ROOT:-mmseqs_chunking}
OUTPUT_ROOT=${OUTPUT_ROOT:-mmseqs}
FASTA_ROOT=${FASTA_ROOT:-fasta}
TARGET_BASES=${TARGET_BASES:-1000000000}
[[ "$TARGET_BASES" =~ ^[1-9][0-9]*$ ]] || { echo 'ERROR: TARGET_BASES must be positive.' >&2; exit 2; }

split_root="$WORK_ROOT/splits"
chunk_root="$WORK_ROOT/chunks"
split_tasks="$WORK_ROOT/split_tasks.tsv"
chunk_tasks="$WORK_ROOT/chunk_tasks.tsv"
merge_tasks="$WORK_ROOT/merge_tasks.tsv"
job_ids="$WORK_ROOT/job_ids.tsv"
mkdir -p "$WORK_ROOT" slurm_output/mmseqs_slurm
tmp_tasks="$WORK_ROOT/.split_tasks.$$.tmp"
trap 'rm -f -- "$tmp_tasks"' EXIT

declare -A samples=()
fasta_end=${FILEEND/.fastq/.fasta}
while IFS= read -r read_file; do
    [[ -n "$read_file" && "$read_file" == *"$FILEEND" ]] || {
        echo "ERROR: invalid reads.txt entry: $read_file" >&2
        exit 2
    }
    sample=${read_file%"$FILEEND"}
    [[ -z ${samples[$sample]+x} ]] || { echo "ERROR: duplicate sample: $sample" >&2; exit 2; }
    samples[$sample]=1
    fasta="$FASTA_ROOT/${read_file%"$FILEEND"}$fasta_end"
    [[ -s "$fasta" ]] || { echo "ERROR: missing FASTA: $fasta" >&2; exit 2; }
    printf '%s\t%s\t%s\n' "$sample" "$fasta" "$TARGET_BASES" >> "$tmp_tasks"
done < reads.txt
sample_count=${#samples[@]}
(( sample_count > 0 )) || { echo 'ERROR: reads.txt contains no samples.' >&2; exit 2; }
mv "$tmp_tasks" "$split_tasks"

split_job=$(sbatch --parsable --array=1-"$sample_count" \
    --export=ALL,MMSEQS_CHUNK_MERGE_DIR="$SCRIPT_DIR",SPLIT_TASKS="$split_tasks",SPLIT_ROOT="$split_root" \
    "$SCRIPT_DIR/split_fasta.slurm")
prepare_job=$(sbatch --parsable --dependency=afterok:"$split_job" \
    --export=ALL,MMSEQS_CHUNK_MERGE_DIR="$SCRIPT_DIR",SPLIT_ROOT="$split_root",CHUNK_ROOT="$chunk_root",OUTPUT_ROOT="$OUTPUT_ROOT",CHUNK_TASKS="$chunk_tasks",MERGE_TASKS="$merge_tasks",JOB_IDS="$job_ids" \
    "$SCRIPT_DIR/prepare_and_submit.slurm")
printf 'split_job=%s\nprepare_job=%s\njob_ids=%s\n' "$split_job" "$prepare_job" "$job_ids"
