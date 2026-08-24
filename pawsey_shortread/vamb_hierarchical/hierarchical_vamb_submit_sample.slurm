#!/bin/bash
# Submit the deterministic hierarchical VAMB pipeline for one sample.
# Run this from the analysis directory containing vamb_crass/SAMPLE/vamb and
# copies or symlinks of the workflow scripts in this directory.

#SBATCH --job-name=hierVAMBsubmit
#SBATCH --partition=work
#SBATCH --account=pawsey1018
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm_output/vamb_slurm/%x-%j.out
#SBATCH --error=slurm_output/vamb_slurm/%x-%j.err

set -euo pipefail

sample=${1:?Usage: hierarchical_vamb_submit_sample.slurm SAMPLE [SHARDS]}
shards=${2:-256}

if ! [[ $shards =~ ^[1-9][0-9]*$ ]]; then
    printf 'SHARDS must be a positive integer, got %q\n' "$shards" >&2
    exit 2
fi

cd "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is required}"

partition_job=$(sbatch --parsable \
    --export="ALL,HIERARCHICAL_STAGE=partition,HIERARCHICAL_SHARDS=$shards" \
    hierarchical_vamb.slurm "$sample")
local_job=$(sbatch --parsable --dependency=afterok:"$partition_job" \
    --array="0-$((shards - 1))" \
    --export="ALL,HIERARCHICAL_STAGE=local,HIERARCHICAL_SHARDS=$shards" \
    hierarchical_vamb.slurm "$sample")
collect_job=$(sbatch --parsable --dependency=afterok:"$local_job" \
    hierarchical_vamb_collect.slurm "$sample")
post_job=$(sbatch --parsable --dependency=afterok:"$collect_job" \
    hierarchical_vamb_submit_reconcile_checkm.slurm "$sample")

printf '%s partition=%s local_array=%s collect=%s post_submit=%s shards=%s\n' \
    "$sample" "$partition_job" "$local_job" "$collect_job" "$post_job" "$shards"
