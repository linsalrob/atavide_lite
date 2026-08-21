#!/bin/bash
# CPU-only representative collection. This stage reads local cluster tables and
# writes the compact merge input; it performs no CUDA/ROCm operation.

#SBATCH --job-name=hierVAMBcollect
#SBATCH --partition=work
#SBATCH --account=pawsey1018
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm_output/vamb_slurm/%x-%j.out
#SBATCH --error=slurm_output/vamb_slurm/%x-%j.err

set -euo pipefail

sample=${1:?Usage: hierarchical_vamb_collect.slurm SAMPLE}
shards=${HIERARCHICAL_SHARDS:-256}
force_args=()
if [[ "${HIERARCHICAL_FORCE_COLLECT:-0}" == 1 ]]; then
    force_args+=(--force)
fi

VENV_DIR="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/vamb"
VAMBDIR="vamb_crass/$sample/vamb"
OUTDIR="${HIERARCHICAL_OUTDIR:-vamb_crass/$sample/hierarchical_vamb}"

[[ -d "$VENV_DIR" ]] || { echo "Missing VAMB environment: $VENV_DIR" >&2; exit 2; }
[[ -f "$VAMBDIR/latent.npz" ]] || { echo "Missing VAMB checkpoint: $VAMBDIR/latent.npz" >&2; exit 3; }

eval "$(conda shell.bash hook)"
conda activate "$VENV_DIR"
python hierarchical_vamb_clustering_rocm.py \
    --vamb-dir "$VAMBDIR" --outdir "$OUTDIR" collect-representatives \
    --shards "$shards" "${force_args[@]}"
