#!/bin/bash
# Submit one CheckM-guided refinement pass after a sample's initial CheckM
# validation has completed. The array size is derived from the QA report.

#SBATCH --job-name=hierVAMBrefSubmit
#SBATCH --partition=work
#SBATCH --account=pawsey1018
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm_output/vamb_slurm/%x-%j.out
#SBATCH --error=slurm_output/vamb_slurm/%x-%j.err

set -euo pipefail

sample=${1:?Usage: hierarchical_vamb_submit_refinement.slurm SAMPLE}
cd "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is required}"
qa="vamb_crass/$sample/hierarchical_vamb/checkm/qa.tsv"
[[ -s "$qa" ]] || { echo "Missing CheckM QA table: $qa" >&2; exit 2; }

count=$(python - "$qa" <<'PY'
import csv
import sys
with open(sys.argv[1], newline="") as handle:
    rows = csv.DictReader(handle, delimiter="\t")
    print(sum(float(row["Completeness"]) >= 90 and float(row["Contamination"]) >= 5 for row in rows))
PY
)

if [[ "$count" -eq 0 ]]; then
    echo "$sample has no bins requiring refinement."
    exit 0
fi

prepare=$(sbatch --parsable hierarchical_vamb_refine_cpu.slurm "$sample" prepare)
split=$(sbatch --parsable --dependency=afterok:"$prepare" --array="0-$((count - 1))" \
    hierarchical_vamb_refine_gpu.slurm "$sample")
merge=$(sbatch --parsable --dependency=afterok:"$split" hierarchical_vamb_refine_cpu.slurm "$sample" merge)
fasta=$(sbatch --parsable --dependency=afterok:"$merge" hierarchical_vamb_refine_cpu.slurm "$sample" fasta)
checkm=$(sbatch --parsable --dependency=afterok:"$fasta" hierarchical_vamb_checkm.slurm \
    "vamb_crass/$sample/hierarchical_vamb/refinement/refined_unsplit_bins" \
    "vamb_crass/$sample/hierarchical_vamb/refinement/checkm")
printf '%s candidates=%s prepare=%s split=%s merge=%s fasta=%s checkm=%s\n' \
    "$sample" "$count" "$prepare" "$split" "$merge" "$fasta" "$checkm"
