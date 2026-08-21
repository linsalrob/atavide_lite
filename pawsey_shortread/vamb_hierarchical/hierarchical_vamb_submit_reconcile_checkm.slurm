#!/bin/bash
# Submit the deterministic post-collection workflow for one sample. This
# lightweight CPU job lets the GPU reconciliation request be made only after
# representative collection has succeeded, then chains export and CheckM.

#SBATCH --job-name=hierVAMBpost
#SBATCH --partition=work
#SBATCH --account=pawsey1018
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm_output/vamb_slurm/%x-%j.out
#SBATCH --error=slurm_output/vamb_slurm/%x-%j.err

set -euo pipefail

sample=${1:?Usage: hierarchical_vamb_submit_reconcile_checkm.slurm SAMPLE}
cd "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is required}"

prepare_job=$(sbatch --parsable --partition=work --account=pawsey1018 \
    hierarchical_vamb_reconcile_cpu.slurm "$sample" prepare)
global_job=$(sbatch --parsable --dependency=afterok:"$prepare_job" \
    hierarchical_vamb_reconcile_gpu.slurm "$sample")
expand_job=$(sbatch --parsable --dependency=afterok:"$global_job" \
    --partition=work --account=pawsey1018 \
    hierarchical_vamb_reconcile_cpu.slurm "$sample" expand)
fasta_job=$(sbatch --parsable --dependency=afterok:"$expand_job" \
    --partition=work --account=pawsey1018 \
    hierarchical_vamb_reconcile_cpu.slurm "$sample" fasta)
checkm_job=$(sbatch --parsable --dependency=afterok:"$fasta_job" \
    hierarchical_vamb_checkm.slurm \
    "vamb_crass/$sample/hierarchical_vamb/hierarchical_unsplit_bins" \
    "vamb_crass/$sample/hierarchical_vamb/checkm")

printf '%s prepare=%s global=%s expand=%s fasta=%s checkm=%s\n' \
    "$sample" "$prepare_job" "$global_job" "$expand_job" "$fasta_job" "$checkm_job"
