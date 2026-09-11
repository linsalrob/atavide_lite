#!/bin/bash
# Materialize medoid-seed children for stable proposals without changing bins.
#SBATCH --job-name=hierVAMBMarkerMat
#SBATCH --partition=work
#SBATCH --account=pawsey1018
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem-per-cpu=1800M
#SBATCH --output=slurm_output/vamb_slurm/%x-%j.out
#SBATCH --error=slurm_output/vamb_slurm/%x-%j.err

# NOTE: Setonix enforces MaxMemPerCPU (work=1840M, highmem=7900M). A bare --mem larger
# than cpus-per-task * MaxMemPerCPU does not fail - Slurm silently raises the CPU count
# to cover it. Request memory with --mem-per-cpu so the two can never disagree.

set -euo pipefail
sample=${1:?Usage: hierarchical_vamb_marker_materialize.slurm SAMPLE}
cd "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is required}"
script_dir=${HIERARCHICAL_SCRIPT_DIR:-"/home/$USER/GitHubs/atavide_lite/pawsey_shortread/vamb_hierarchical"}
VENV_DIR="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/vamb"
root="vamb_crass/$sample/hierarchical_vamb"
eval "$(conda shell.bash hook)"
conda activate "$VENV_DIR"
python "$script_dir/validate_hierarchical_vamb_proposals.py" materialize \
    --decisions "$root/gated_refinement_evaluation/split_decisions.tsv" \
    --parent-bins "$root/hierarchical_unsplit_bins" \
    --output "$root/marker_validation"
