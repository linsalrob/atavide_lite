#!/bin/bash
# Task file: label<TAB>input_fasta<TAB>output_directory

#SBATCH --job-name=mmseqs_chunk
#SBATCH --partition=work
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=1800M
#SBATCH --output=slurm_output/mmseqs_slurm/chunk-%A_%a.out
#SBATCH --error=slurm_output/mmseqs_slurm/chunk-%A_%a.err

# NOTE: Setonix enforces MaxMemPerCPU (work=1840M, highmem=7900M). A bare --mem larger
# than cpus-per-task * MaxMemPerCPU does not fail - Slurm silently raises the CPU count
# to cover it. Request memory with --mem-per-cpu so the two can never disagree.

set -euo pipefail
TASKS=${CHUNK_TASKS:?Set CHUNK_TASKS to a tab-delimited task file.}
task_index=${SLURM_ARRAY_TASK_ID:?Submit as an array.}
line=$(sed -n "${task_index}p" "$TASKS")
IFS=$'\t' read -r label input output <<< "$line"
[[ -n "$label" && -s "$input" && -n "$output" ]] || {
    echo 'ERROR: invalid chunk task.' >&2
    exit 2
}
if [[ -s "$output/DONE" ]]; then
    for suffix in lca.tsv report tophit_aln tophit_report.gz; do
        [[ -s "$output/${label}_${suffix}" ]] || {
            echo "ERROR: $output/DONE exists but ${label}_${suffix} is missing" >&2
            exit 2
        }
    done
    echo "$label already complete"
    exit 0
fi
[[ ! -e "$output" ]] || {
    echo "ERROR: incomplete chunk output exists; inspect before retrying: $output" >&2
    exit 2
}

MMSEQS=${MMSEQS_BIN:-/scratch/${PAWSEY_PROJECT}/${USER}/software/miniconda3/atavide_lite/bin/mmseqs}
TARGET=${TARGET_DB:-/scratch/${PAWSEY_PROJECT}/${USER}/Databases/UniRef100/UniRef100}
[[ -x "$MMSEQS" && -s "${TARGET}.dbtype" ]] || {
    echo 'ERROR: MMseqs or target database is unavailable.' >&2
    exit 2
}

mkdir -p "$output"
tmp="$output/tmp"
"$MMSEQS" easy-taxonomy "$input" "$TARGET" "$output/$label" "$tmp" \
    --threads "${SLURM_CPUS_PER_TASK:-1}" \
    --split-memory-limit "${MMSEQS_SPLIT_MEMORY_LIMIT:-48G}" \
    --remove-tmp-files 0

mapfile -t swapped < <(find "$tmp" -type f -name result_aln_swapped.dbtype -printf '%h/result_aln_swapped\n')
mapfile -t result < <(find "$tmp" -type f -name result.dbtype -printf '%h/result\n')
(( ${#swapped[@]} == 1 && ${#result[@]} == 1 )) || {
    echo "ERROR: expected one retained alignment and taxonomy database for $label" >&2
    exit 2
}
for prefix in "${swapped[0]}" "${result[0]}"; do
    [[ -s "${prefix}.dbtype" && -s "${prefix}.index" ]] || {
        echo "ERROR: incomplete retained MMseqs database: $prefix" >&2
        exit 2
    }
done
for suffix in lca.tsv report tophit_aln tophit_report; do
    [[ -s "$output/${label}_${suffix}" ]] || {
        echo "ERROR: missing easy-taxonomy output: ${label}_${suffix}" >&2
        exit 2
    }
done
pigz -p "${SLURM_CPUS_PER_TASK:-1}" "$output/${label}_tophit_report"
printf 'label=%s\ninput=%s\nresult=%s\nresult_aln_swapped=%s\n' \
    "$label" "$input" "${result[0]}" "${swapped[0]}" > "$output/DONE"
