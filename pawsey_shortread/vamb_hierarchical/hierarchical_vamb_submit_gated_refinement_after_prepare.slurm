#!/bin/bash
# Submit GPU proposal tasks only after their manifest has been prepared.
#SBATCH --job-name=hierVAMBGateContinue
#SBATCH --partition=work
#SBATCH --account=pawsey1018
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm_output/vamb_slurm/%x-%j.out
#SBATCH --error=slurm_output/vamb_slurm/%x-%j.err
set -euo pipefail
sample=${1:?Usage: hierarchical_vamb_submit_gated_refinement_after_prepare.slurm SAMPLE}
cd "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is required}"
script_dir=${HIERARCHICAL_SCRIPT_DIR:-"/home/$USER/GitHubs/atavide_lite/pawsey_shortread/vamb_hierarchical"}
tasks="vamb_crass/$sample/hierarchical_vamb/gated_refinement_evaluation/tasks.tsv"
[[ -s "$tasks" ]] || { echo "Missing proposal task manifest: $tasks" >&2; exit 2; }
count=$(tail -n +2 "$tasks" | wc -l)
if [[ "$count" -eq 0 ]]; then
  echo "$sample has no bounded split proposals; inspect MEGABIN_REVIEW records."
  exit 0
fi
split=$(sbatch --parsable --array="0-$((count - 1))" "$script_dir/hierarchical_vamb_gated_refine_gpu.slurm" "$sample")
summary=$(sbatch --parsable --dependency=afterok:"$split" "$script_dir/hierarchical_vamb_gated_refine_cpu.slurm" "$sample" summarize)
printf '%s proposals=%s split=%s summarize=%s\n' "$sample" "$count" "$split" "$summary"
