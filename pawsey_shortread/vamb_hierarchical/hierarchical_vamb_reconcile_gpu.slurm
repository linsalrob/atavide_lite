#!/bin/bash
# GPU-only global clustering of multi-contig local-bin representatives.
#SBATCH --job-name=hierVAMBreconGPU
#SBATCH --partition=gpu
#SBATCH --account=pawsey1018-gpu
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --output=slurm_output/vamb_slurm/%x-%j.out
#SBATCH --error=slurm_output/vamb_slurm/%x-%j.err
set -euo pipefail
sample=${1:?Usage: hierarchical_vamb_reconcile_gpu.slurm SAMPLE}
VENV_DIR="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/vamb"
VAMBDIR="vamb_crass/$sample"
OUTDIR="$VAMBDIR/hierarchical_vamb"
module load rocm/6.2.4
eval "$(conda shell.bash hook)"
conda activate "$VENV_DIR"
python hierarchical_vamb_clustering_rocm.py --vamb-dir "$VAMBDIR/vamb" --outdir "$OUTDIR" cluster-global
