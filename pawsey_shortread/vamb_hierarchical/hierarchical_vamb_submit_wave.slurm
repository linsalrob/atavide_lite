#!/bin/bash
# Submit later hierarchical-VAMB waves after the preceding wave has fully
# collected its representatives. This runs on the CPU work queue so it does
# not consume a GPU or count against the GPU association's submit limit.
#
# The first job is submitted with a dependency on the initial three samples:
#   sbatch --dependency=afterok:47380672:47380675:47380678 \
#     hierarchical_vamb_submit_wave.slurm 0

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

wave=${1:?Usage: hierarchical_vamb_submit_wave.slurm WAVE_NUMBER}
cd "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is required}"

case "$wave" in
    0) samples=(Eagle ExmouthGulfWater Fiddler) ;;
    1) samples=(Guitarfish NingalooWater Port) ;;
    2) samples=(Seven TantabiddiWater Water) ;;
    3) samples=(Whitespotted) ;;
    *) echo "All deferred hierarchical VAMB waves have been submitted."; exit 0 ;;
esac

collect_jobs=()
for sample in "${samples[@]}"; do
    partition_job=$(sbatch --parsable --partition=gpu --gres=gpu:1 --time=04:00:00 \
        --export=ALL,HIERARCHICAL_STAGE=partition,HIERARCHICAL_SHARDS=256 \
        hierarchical_vamb.slurm "$sample")
    local_job=$(sbatch --parsable --dependency=afterok:"$partition_job" --array=0-255 \
        --partition=gpu --gres=gpu:1 --time=1-00:00:00 \
        --export=ALL,HIERARCHICAL_STAGE=local,HIERARCHICAL_SHARDS=256 \
        hierarchical_vamb.slurm "$sample")
    collect_job=$(sbatch --parsable --dependency=afterok:"$local_job" \
        --export=ALL,HIERARCHICAL_SHARDS=256 \
        hierarchical_vamb_collect.slurm "$sample")
    reconcile_job=$(sbatch --parsable --dependency=afterok:"$collect_job" \
        hierarchical_vamb_submit_reconcile_checkm.slurm "$sample")
    collect_jobs+=("$collect_job")
    printf '%s partition=%s local=%s collect=%s reconcile_submit=%s\n' \
        "$sample" "$partition_job" "$local_job" "$collect_job" "$reconcile_job"
done

if [[ "$wave" -lt 3 ]]; then
    dependencies=$(IFS=:; echo "${collect_jobs[*]}")
    next_wave=$((wave + 1))
    next_job=$(sbatch --parsable --dependency=afterok:"$dependencies" \
        hierarchical_vamb_submit_wave.slurm "$next_wave")
    printf 'next_wave=%s job=%s dependency=afterok:%s\n' "$next_wave" "$next_job" "$dependencies"
fi
