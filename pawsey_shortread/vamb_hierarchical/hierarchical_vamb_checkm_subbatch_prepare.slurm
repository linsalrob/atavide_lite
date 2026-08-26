#!/bin/bash
# Extract one failed parent CheckM batch into a new manifest for smaller retries.
#SBATCH --job-name=checkmSubPrep
#SBATCH --partition=work
#SBATCH --account=pawsey1018
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=slurm_output/checkm-subbatch-prep-%j.out
#SBATCH --error=slurm_output/checkm-subbatch-prep-%j.err

set -euo pipefail
parent_report=${1:?Usage: hierarchical_vamb_checkm_subbatch_prepare.slurm PARENT_REPORT PARENT_INDEX PARENT_BATCH_SIZE SUB_REPORT}
parent_index=${2:?Missing parent batch index}
parent_batch_size=${3:?Missing parent batch size}
sub_report=${4:?Missing subbatch report directory}
parent_manifest="${parent_report}_batches/bins.txt"
sub_manifest_dir="${sub_report}_batches"
sub_manifest="$sub_manifest_dir/bins.txt"
[[ -s "$parent_manifest" ]] || { echo "Missing parent manifest: $parent_manifest" >&2; exit 2; }
[[ "$parent_index" =~ ^[0-9]+$ ]] || { echo "Invalid parent index: $parent_index" >&2; exit 2; }
[[ "$parent_batch_size" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid parent batch size: $parent_batch_size" >&2; exit 2; }
[[ ! -e "$sub_manifest" ]] || { echo "Refusing to overwrite $sub_manifest" >&2; exit 3; }
mkdir -p "$sub_manifest_dir" "$sub_report"
start=$((parent_index * parent_batch_size + 1))
end=$((start + parent_batch_size - 1))
sed -n "${start},${end}p" "$parent_manifest" > "$sub_manifest"
count=$(wc -l < "$sub_manifest")
[[ "$count" -gt 0 ]] || { echo "No manifest rows in parent batch $parent_index" >&2; exit 4; }
printf 'parent_batch=%s rows=%s-%s bins=%s manifest=%s\n' "$parent_index" "$start" "$end" "$count" "$sub_manifest"
