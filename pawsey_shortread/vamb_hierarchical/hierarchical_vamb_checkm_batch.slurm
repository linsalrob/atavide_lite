#!/bin/bash
# Validate one independent batch of bins using node-local CheckM work space.
#SBATCH --job-name=checkmBatch
#SBATCH --partition=highmem
#SBATCH --time=2-0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB
#SBATCH --output=slurm_output/checkm-batch-%A_%a.out
#SBATCH --error=slurm_output/checkm-batch-%A_%a.err
set -euo pipefail
source_dir=${1:?Usage: hierarchical_vamb_checkm_batch.slurm SOURCE_BINS REPORT_DIR BATCH_SIZE}
report_dir=${2:?Usage: hierarchical_vamb_checkm_batch.slurm SOURCE_BINS REPORT_DIR BATCH_SIZE}
batch_size=${3:?Usage: hierarchical_vamb_checkm_batch.slurm SOURCE_BINS REPORT_DIR BATCH_SIZE}
: "${SLURM_ARRAY_TASK_ID:?Requires a Slurm array task id}"
manifest_dir="${report_dir}_batches"
manifest="$manifest_dir/bins.txt"
start=$((SLURM_ARRAY_TASK_ID * batch_size + 1))
end=$((start + batch_size - 1))
work_dir="${SLURM_TMPDIR:-/tmp}/checkm-${SLURM_JOB_ID}"
input_dir="$work_dir/input"
output_dir="$work_dir/output"
mkdir -p "$input_dir" "$output_dir/tmp"
sed -n "${start},${end}p" "$manifest" | while IFS= read -r fasta; do
    ln -s "$fasta" "$input_dir/$(basename "$fasta")"
done
[[ -n $(find "$input_dir" -maxdepth 1 -type l -print -quit) ]] || {
    echo "No bins assigned to task $SLURM_ARRAY_TASK_ID" >&2
    exit 2
}

eval "$(conda shell.bash hook)"
"$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh"
ATAVIDE_CONDA="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite"
conda activate "$ATAVIDE_CONDA"
export CHECKM_DATA_PATH="${CHECKM_DATA_PATH:-$ATAVIDE_CONDA/checkm_data}"
site_dir="$work_dir/sitecustomize"
mkdir -p "$site_dir"
printf '%s\n' 'import sys' 'sys.setrecursionlimit(max(sys.getrecursionlimit(), 50000))' > "$site_dir/sitecustomize.py"
export PYTHONPATH="$site_dir${PYTHONPATH:+:$PYTHONPATH}"

checkm lineage_wf --tmpdir "$output_dir/tmp" -t 16 "$input_dir" "$output_dir"
checkm qa --tab_table -o 1 -f "$manifest_dir/qa_${SLURM_ARRAY_TASK_ID}.tsv" \
    --tmpdir "$output_dir/tmp" -t 16 "$output_dir/lineage.ms" "$output_dir"
