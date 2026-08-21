#!/bin/bash
# Stand-alone hierarchical VAMB workflow. Inputs are a completed VAMB feature
# directory with composition.npz, abundance.npz, and latent.npz; no prior
# cluster table or resume_vamb_clustering_rocm.py output is required.
#
# Stage 1 (one GPU; use gpu-dev for a pilot):
#   sbatch --partition=gpu-dev --gres=gpu:1 --time=4:00:00 \
#     --export=ALL,HIERARCHICAL_STAGE=partition,HIERARCHICAL_SHARDS=256 \
#     hierarchical_vamb.slurm Whaleshark
#
# Stage 2 (one GPU per array task after stage 1 succeeds):
#   sbatch --array=0-255%16 --gres=gpu:1 --time=1-00:00:00 \
#     --export=ALL,HIERARCHICAL_STAGE=local,HIERARCHICAL_SHARDS=256 \
#     hierarchical_vamb.slurm Whaleshark
#
# Stage 3 (concatenate representatives for global reconciliation):
#   sbatch --export=ALL,HIERARCHICAL_STAGE=collect,HIERARCHICAL_SHARDS=256 \
#     hierarchical_vamb.slurm Whaleshark

#SBATCH --job-name=hierVAMB
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=gpu
#SBATCH --account=pawsey1018-gpu
#SBATCH --gres=gpu:1
#SBATCH -o slurm_output/vamb_slurm/%x-%A_%a.out
#SBATCH -e slurm_output/vamb_slurm/%x-%A_%a.err

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: sbatch ... hierarchical_vamb.slurm SAMPLE" >&2
    exit 2
fi

SAMPLE=$1
STAGE=${HIERARCHICAL_STAGE:-partition}
SHARDS=${HIERARCHICAL_SHARDS:-256}
VAMBDIR="vamb_crass/$SAMPLE/vamb"
OUTDIR="${HIERARCHICAL_OUTDIR:-vamb_crass/$SAMPLE/hierarchical_vamb}"
VENV_DIR="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/vamb"

for required in "$VAMBDIR/composition.npz" "$VAMBDIR/abundance.npz" "$VAMBDIR/latent.npz"; do
    [[ -f "$required" ]] || { echo "Missing required VAMB feature input: $required" >&2; exit 3; }
done
[[ -d "$VENV_DIR" ]] || { echo "Missing VAMB environment: $VENV_DIR" >&2; exit 4; }

module load rocm/6.2.4
eval "$(conda shell.bash hook)"
conda activate "$VENV_DIR"

echo "[$(date -Is)] stage=$STAGE sample=$SAMPLE shards=$SHARDS outdir=$OUTDIR" >&2

case "$STAGE" in
    partition)
        python hierarchical_vamb_clustering_rocm.py \
            --vamb-dir "$VAMBDIR" --outdir "$OUTDIR" partition \
            --shards "$SHARDS"
        ;;
    local)
        : "${SLURM_ARRAY_TASK_ID:?local stage requires an array task ID}"
        shard_dir="$OUTDIR/$(printf 'shard_%04d' "$SLURM_ARRAY_TASK_ID")"
        if [[ -f "$shard_dir/local_clusters.tsv" && -f "$shard_dir/representatives.npz" ]]; then
            echo "[$(date -Is)] shard $SLURM_ARRAY_TASK_ID already complete; skipping" >&2
            exit 0
        fi
        python hierarchical_vamb_clustering_rocm.py \
            --vamb-dir "$VAMBDIR" --outdir "$OUTDIR" local \
            --shard "$SLURM_ARRAY_TASK_ID"
        ;;
    collect)
        collect_args=()
        if [[ "${HIERARCHICAL_FORCE_COLLECT:-0}" == 1 ]]; then
            collect_args+=(--force)
        fi
        python hierarchical_vamb_clustering_rocm.py \
            --vamb-dir "$VAMBDIR" --outdir "$OUTDIR" collect-representatives \
            --shards "$SHARDS" "${collect_args[@]}"
        ;;
    *)
        echo "Unknown HIERARCHICAL_STAGE: $STAGE (partition, local, or collect)" >&2
        exit 5
        ;;
esac
