#!/bin/bash
# Prepare and summarize non-destructive, multi-seed split proposals.
#SBATCH --job-name=hierVAMBGateCPU
#SBATCH --partition=work
#SBATCH --account=pawsey1018
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --output=slurm_output/vamb_slurm/%x-%j.out
#SBATCH --error=slurm_output/vamb_slurm/%x-%j.err
set -euo pipefail
sample=${1:?Usage: hierarchical_vamb_gated_refine_cpu.slurm SAMPLE [prepare|summarize]}
stage=${2:?Specify prepare or summarize}
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
VENV_DIR="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/vamb"
sample_dir="vamb_crass/$sample"
outdir="$sample_dir/hierarchical_vamb"
qa_path="${EVALUATION_QA:-$outdir/checkm/qa.tsv}"
eval "$(conda shell.bash hook)"
conda activate "$VENV_DIR"
case "$stage" in
  prepare) python "$script_dir/evaluate_hierarchical_vamb_refinement.py" --vamb-dir "$sample_dir/vamb" --outdir "$outdir" prepare --qa "$qa_path" --assignments "$outdir/hierarchical_unsplit.tsv" ;;
  summarize) python "$script_dir/evaluate_hierarchical_vamb_refinement.py" --vamb-dir "$sample_dir/vamb" --outdir "$outdir" summarize ;;
  *) echo "Unknown stage: $stage" >&2; exit 2 ;;
esac
