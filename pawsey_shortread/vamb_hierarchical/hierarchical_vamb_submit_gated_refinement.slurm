#!/bin/bash
# Submit proposal generation and stability assessment only; no automatic bin replacement.
#SBATCH --job-name=hierVAMBGateSubmit
#SBATCH --partition=work
#SBATCH --account=pawsey1018
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm_output/vamb_slurm/%x-%j.out
#SBATCH --error=slurm_output/vamb_slurm/%x-%j.err
set -euo pipefail
sample=${1:?Usage: hierarchical_vamb_submit_gated_refinement.slurm SAMPLE}
cd "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is required}"
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
prepare=$(sbatch --parsable "$script_dir/hierarchical_vamb_gated_refine_cpu.slurm" "$sample" prepare)
continue=$(sbatch --parsable --dependency=afterok:"$prepare" "$script_dir/hierarchical_vamb_submit_gated_refinement_after_prepare.slurm" "$sample")
printf '%s prepare=%s continue=%s\n' "$sample" "$prepare" "$continue"
