#!/bin/bash

#SBATCH --job-name=prepare_mmseqs
#SBATCH --partition=work
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=slurm_output/mmseqs_slurm/prepare-%j.out
#SBATCH --error=slurm_output/mmseqs_slurm/prepare-%j.err

set -euo pipefail
: "${SPLIT_ROOT:?Set SPLIT_ROOT.}"
: "${CHUNK_ROOT:?Set CHUNK_ROOT.}"
: "${OUTPUT_ROOT:?Set OUTPUT_ROOT.}"
: "${CHUNK_TASKS:?Set CHUNK_TASKS.}"
: "${MERGE_TASKS:?Set MERGE_TASKS.}"
: "${JOB_IDS:?Set JOB_IDS.}"
SCRIPT_DIR=${MMSEQS_CHUNK_MERGE_DIR:-$HOME/GitHubs/atavide_lite/pawsey_minion/mmseqs_chunk_merge}
[[ -s "$SCRIPT_DIR/build_tasks.py" ]] || {
    echo "ERROR: MMseqs chunk/merge module unavailable: $SCRIPT_DIR" >&2
    exit 2
}

/usr/bin/python3.11 "$SCRIPT_DIR/build_tasks.py" \
    --split-root "$SPLIT_ROOT" --chunk-root "$CHUNK_ROOT" --output-root "$OUTPUT_ROOT" \
    --chunk-tasks "$CHUNK_TASKS" --merge-tasks "$MERGE_TASKS"
chunk_count=$(wc -l < "$CHUNK_TASKS")
merge_count=$(wc -l < "$MERGE_TASKS")
(( chunk_count > 0 && merge_count > 0 )) || { echo 'ERROR: empty task manifest.' >&2; exit 2; }

chunk_job=$(sbatch --parsable --array=1-"$chunk_count" \
    --export=ALL,CHUNK_TASKS="$CHUNK_TASKS" "$SCRIPT_DIR/search_chunk.slurm")

declare -A dependencies=()
index=0
while IFS=$'\t' read -r label _; do
    sample=${label%%__*}
    ((index += 1))
    dependencies["$sample"]+="${chunk_job}_${index}:"
done < "$CHUNK_TASKS"

merge_jobs=()
index=0
while IFS=$'\t' read -r sample _; do
    ((index += 1))
    dependency=${dependencies[$sample]%:}
    merge_job=$(sbatch --parsable --array="$index-$index" --dependency=afterok:"$dependency" \
        --export=ALL,MERGE_TASKS="$MERGE_TASKS" "$SCRIPT_DIR/merge_sample.slurm")
    merge_jobs+=("$merge_job")
done < "$MERGE_TASKS"

merge_dependency=$(IFS=:; echo "${merge_jobs[*]}")
{
    printf 'chunk_job=%s\nchunk_count=%s\n' "$chunk_job" "$chunk_count"
    printf 'merge_count=%s\nmerge_jobs=%s\n' "$merge_count" "$merge_dependency"
} > "$JOB_IDS"
printf 'Submitted %s chunks as %s and %s sample merges: %s\n' \
    "$chunk_count" "$chunk_job" "$merge_count" "$merge_dependency"
