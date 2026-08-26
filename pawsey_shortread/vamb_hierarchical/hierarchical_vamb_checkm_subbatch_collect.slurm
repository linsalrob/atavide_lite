#!/bin/bash
# Merge successful subbatches into the missing QA table for one parent batch.
#SBATCH --job-name=checkmSubCollect
#SBATCH --partition=work
#SBATCH --account=pawsey1018
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=slurm_output/checkm-subbatch-collect-%j.out
#SBATCH --error=slurm_output/checkm-subbatch-collect-%j.err

set -euo pipefail
sub_report=${1:?Usage: hierarchical_vamb_checkm_subbatch_collect.slurm SUB_REPORT SUBBATCH_COUNT PARENT_REPORT PARENT_INDEX}
subbatch_count=${2:?Missing subbatch count}
parent_report=${3:?Missing parent report directory}
parent_index=${4:?Missing parent batch index}
sub_manifest_dir="${sub_report}_batches"
parent_table="${parent_report}_batches/qa_${parent_index}.tsv"
temporary="${parent_table}.tmp-${SLURM_JOB_ID}"
[[ "$subbatch_count" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid subbatch count: $subbatch_count" >&2; exit 2; }
[[ "$parent_index" =~ ^[0-9]+$ ]] || { echo "Invalid parent index: $parent_index" >&2; exit 2; }
[[ ! -e "$parent_table" ]] || { echo "Refusing to overwrite $parent_table" >&2; exit 3; }
first="$sub_manifest_dir/qa_0.tsv"
[[ -s "$first" ]] || { echo "Missing $first" >&2; exit 4; }
trap 'rm -f "$temporary"' EXIT
head -n 1 "$first" > "$temporary"
for task in $(seq 0 "$((subbatch_count - 1))"); do
    table="$sub_manifest_dir/qa_${task}.tsv"
    [[ -s "$table" ]] || { echo "Missing $table" >&2; exit 5; }
    [[ "$(head -n 1 "$table")" == "$(head -n 1 "$first")" ]] || { echo "Header mismatch in $table" >&2; exit 6; }
    tail -n +2 "$table" >> "$temporary"
done
mv "$temporary" "$parent_table"
trap - EXIT
printf 'Merged %s subbatches into %s (%s assessed bins)\n' "$subbatch_count" "$parent_table" "$(( $(wc -l < "$parent_table") - 1 ))"
