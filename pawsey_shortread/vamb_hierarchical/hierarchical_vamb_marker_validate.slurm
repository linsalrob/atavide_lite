#!/bin/bash
# Re-establish each parent lineage marker set, then apply that same set to its
# candidate children and retain compact marker IDs/copy locations.
#SBATCH --job-name=hierVAMBMarkerQC
#SBATCH --partition=highmem
#SBATCH --account=pawsey1018
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --output=slurm_output/vamb_slurm/%x-%j.out
#SBATCH --error=slurm_output/vamb_slurm/%x-%j.err

set -euo pipefail
sample=${1:?Usage: hierarchical_vamb_marker_validate.slurm SAMPLE}
cd "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is required}"
script_dir=${HIERARCHICAL_SCRIPT_DIR:-"/home/$USER/GitHubs/atavide_lite/pawsey_shortread/vamb_hierarchical"}
root="vamb_crass/$sample/hierarchical_vamb/marker_validation"
report="$root/checkm_parent_matched"
work="${SLURM_TMPDIR:-/tmp}/marker-validation-${SLURM_JOB_ID}"
parent_work="$work/parents"
matched_work="$work/parent-matched"
ATAVIDE_CONDA="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite"

[[ -s "$root/bin_manifest.tsv" ]] || { echo "Missing marker-validation manifest" >&2; exit 2; }
if [[ -e "$root/NO_STABLE_PROPOSALS" ]]; then
    echo "$sample has no proposals passing both stability gates; retaining all parents."
    exit 0
fi
[[ ! -e "$report" ]] || { echo "Refusing to overwrite $report" >&2; exit 3; }
mkdir -p "$work/tmp" "$report"
eval "$(conda shell.bash hook)"
"$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh"
conda activate "$ATAVIDE_CONDA"
export CHECKM_DATA_PATH="${CHECKM_DATA_PATH:-$ATAVIDE_CONDA/checkm_data}"
export MPLCONFIGDIR="$work/matplotlib"
site_dir="$work/sitecustomize"
mkdir -p "$site_dir" "$MPLCONFIGDIR"
printf '%s\n' 'import sys' 'sys.setrecursionlimit(max(sys.getrecursionlimit(), 50000))' > "$site_dir/sitecustomize.py"
export PYTHONPATH="$site_dir${PYTHONPATH:+:$PYTHONPATH}"

checkm lineage_wf --tmpdir "$work/tmp" -t 16 "$root/parents" "$parent_work"
cp "$parent_work/lineage.ms" "$report/parent_lineage.ms"

conda activate "/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/vamb"
python "$script_dir/validate_hierarchical_vamb_proposals.py" build-marker-file \
    --parent-lineage "$report/parent_lineage.ms" \
    --manifest "$root/bin_manifest.tsv" \
    --output "$report/parent_matched.ms"

conda activate "$ATAVIDE_CONDA"
checkm analyze --tmpdir "$work/tmp" -t 16 "$report/parent_matched.ms" "$root/bins" "$matched_work"
checkm qa --tab_table -o 1 -f "$report/qa.tsv" --tmpdir "$work/tmp" -t 16 \
    "$report/parent_matched.ms" "$matched_work"
cp "$matched_work/storage/marker_gene_stats.tsv" "$report/marker_gene_stats.tsv"
cp "$matched_work/storage/bin_stats_ext.tsv" "$report/bin_stats_ext.tsv"

conda activate "/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/vamb"
python "$script_dir/validate_hierarchical_vamb_proposals.py" report \
    --manifest "$root/bin_manifest.tsv" \
    --marker-stats "$report/marker_gene_stats.tsv" \
    --qa "$report/qa.tsv" \
    --conflicts "$root/marker_conflicts.tsv" \
    --output "$root/marker_gate_decisions.tsv"
