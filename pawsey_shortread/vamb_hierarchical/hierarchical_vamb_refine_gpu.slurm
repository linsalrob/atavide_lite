#!/bin/bash
# One CheckM-selected parent bin per GPU array task.
#SBATCH --job-name=hierVAMBrefGPU
#SBATCH --partition=gpu
#SBATCH --account=pawsey1018-gpu
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --output=slurm_output/vamb_slurm/%x-%A_%a.out
#SBATCH --error=slurm_output/vamb_slurm/%x-%A_%a.err
set -euo pipefail
sample=${1:?Usage: hierarchical_vamb_refine_gpu.slurm SAMPLE}
: "${SLURM_ARRAY_TASK_ID:?Requires a Slurm array task id}"
VENV_DIR="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/vamb"
SAMPLE_DIR="vamb_crass/$sample"
OUTDIR="$SAMPLE_DIR/hierarchical_vamb"
module load rocm/6.2.4
eval "$(conda shell.bash hook)"
conda activate "$VENV_DIR"
python refine_hierarchical_vamb.py --vamb-dir "$SAMPLE_DIR/vamb" --outdir "$OUTDIR" split --task-id "$SLURM_ARRAY_TASK_ID"
