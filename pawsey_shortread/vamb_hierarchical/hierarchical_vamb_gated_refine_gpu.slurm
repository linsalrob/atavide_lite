#!/bin/bash
# One k/seed split proposal per GPU array task; no bin assignments are replaced.
#SBATCH --job-name=hierVAMBGateGPU
#SBATCH --partition=gpu
#SBATCH --account=pawsey1018-gpu
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --output=slurm_output/vamb_slurm/%x-%A_%a.out
#SBATCH --error=slurm_output/vamb_slurm/%x-%A_%a.err
set -euo pipefail
sample=${1:?Usage: hierarchical_vamb_gated_refine_gpu.slurm SAMPLE}
: "${SLURM_ARRAY_TASK_ID:?Requires an array task id}"
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
VENV_DIR="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/vamb"
sample_dir="vamb_crass/$sample"
module load rocm/6.2.4
eval "$(conda shell.bash hook)"
conda activate "$VENV_DIR"
python "$script_dir/evaluate_hierarchical_vamb_refinement.py" --vamb-dir "$sample_dir/vamb" --outdir "$sample_dir/hierarchical_vamb" split --task-id "$SLURM_ARRAY_TASK_ID"
