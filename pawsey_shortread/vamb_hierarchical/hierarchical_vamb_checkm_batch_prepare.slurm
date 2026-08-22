#!/bin/bash
# Prepare a stable, lightweight manifest for batched CheckM validation.
#SBATCH --job-name=checkmBatchPrep
#SBATCH --partition=work
#SBATCH --account=pawsey1018
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=slurm_output/checkm-batch-prep-%j.out
#SBATCH --error=slurm_output/checkm-batch-prep-%j.err
set -euo pipefail
source_dir=${1:?Usage: hierarchical_vamb_checkm_batch_prepare.slurm SOURCE_BINS REPORT_DIR}
report_dir=${2:?Usage: hierarchical_vamb_checkm_batch_prepare.slurm SOURCE_BINS REPORT_DIR}
manifest_dir="${report_dir}_batches"
mkdir -p "$manifest_dir" "$report_dir"
find "$source_dir" -maxdepth 1 -type f -name '*.fna' -printf '%p\n' | sort > "$manifest_dir/bins.txt"
count=$(wc -l < "$manifest_dir/bins.txt")
[[ "$count" -gt 0 ]] || { echo "No FASTA bins in $source_dir" >&2; exit 2; }
printf 'bins=%s\nmanifest=%s\n' "$count" "$manifest_dir/bins.txt"
