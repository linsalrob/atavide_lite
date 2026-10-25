#!/bin/bash
# Task file: sample<TAB>input_fasta<TAB>target_bases

#SBATCH --job-name=mmseqs_split
#SBATCH --partition=work
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=1800M
#SBATCH --output=slurm_output/mmseqs_slurm/split-%A_%a.out
#SBATCH --error=slurm_output/mmseqs_slurm/split-%A_%a.err

# NOTE: Setonix enforces MaxMemPerCPU (work=1840M, highmem=7900M). A bare --mem larger
# than cpus-per-task * MaxMemPerCPU does not fail - Slurm silently raises the CPU count
# to cover it. Request memory with --mem-per-cpu so the two can never disagree.

set -euo pipefail
TASKS=${SPLIT_TASKS:?Set SPLIT_TASKS to a tab-delimited task file.}
task_index=${SLURM_ARRAY_TASK_ID:?Submit as an array.}
line=$(sed -n "${task_index}p" "$TASKS")
IFS=$'\t' read -r sample input target_bases <<< "$line"
[[ -n "$sample" && -s "$input" && "$target_bases" =~ ^[1-9][0-9]*$ ]] || {
    echo 'ERROR: invalid split task.' >&2
    exit 2
}

SPLIT_ROOT=${SPLIT_ROOT:-mmseqs_chunking/splits}
output="$SPLIT_ROOT/$sample"
[[ ! -s "$output/DONE" ]] || { echo "$sample split already complete"; exit 0; }
[[ ! -e "$output" ]] || {
    echo "ERROR: incomplete split directory exists; inspect before retrying: $output" >&2
    exit 2
}

SCRIPT_DIR=${MMSEQS_CHUNK_MERGE_DIR:-$HOME/GitHubs/atavide_lite/pawsey_minion/mmseqs_chunk_merge}
[[ -s "$SCRIPT_DIR/split_fasta_by_bases.py" ]] || {
    echo "ERROR: MMseqs chunk/merge module unavailable: $SCRIPT_DIR" >&2
    exit 2
}
python3 "$SCRIPT_DIR/split_fasta_by_bases.py" \
    --input "$input" --output-dir "$output" --target-bases "$target_bases"
