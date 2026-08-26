#!/bin/bash
# Merge batched CheckM QA tables into the standard downstream location.
#SBATCH --job-name=checkmBatchCollect
#SBATCH --partition=work
#SBATCH --account=pawsey1018
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=slurm_output/checkm-batch-collect-%j.out
#SBATCH --error=slurm_output/checkm-batch-collect-%j.err
set -euo pipefail
report_dir=${1:?Usage: hierarchical_vamb_checkm_batch_collect.slurm REPORT_DIR BATCH_COUNT}
batch_count=${2:?Usage: hierarchical_vamb_checkm_batch_collect.slurm REPORT_DIR BATCH_COUNT}
manifest_dir="${report_dir}_batches"
output="$report_dir/qa.tsv"
[[ ! -e "$output" ]] || { echo "Refusing to overwrite $output" >&2; exit 2; }
first="$manifest_dir/qa_0.tsv"
[[ -s "$first" ]] || { echo "Missing $first" >&2; exit 3; }
head -n 1 "$first" > "$output"
for task in $(seq 0 "$((batch_count - 1))"); do
    table="$manifest_dir/qa_${task}.tsv"
    [[ -s "$table" ]] || { echo "Missing $table" >&2; exit 4; }
    tail -n +2 "$table" >> "$output"
done
printf 'Merged %s QA batches into %s (%s assessed bins)\n' "$batch_count" "$output" "$(( $(wc -l < "$output") - 1 ))"
