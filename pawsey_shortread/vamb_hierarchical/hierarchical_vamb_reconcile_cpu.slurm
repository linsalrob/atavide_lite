#!/bin/bash
# CPU stages of global reconciliation: prepare candidates, expand memberships,
# and write FASTA bins after expansion.
#SBATCH --job-name=hierVAMBreconCPU
#SBATCH --partition=work
#SBATCH --account=pawsey1018
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=slurm_output/vamb_slurm/%x-%j.out
#SBATCH --error=slurm_output/vamb_slurm/%x-%j.err
set -euo pipefail
sample=${1:?Usage: hierarchical_vamb_reconcile_cpu.slurm SAMPLE [prepare|expand|fasta]}
stage=${2:?Specify prepare, expand, or fasta}
VENV_DIR="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/vamb"
SAMPLE_DIR="vamb_crass/$sample"
OUTDIR="$SAMPLE_DIR/hierarchical_vamb"
eval "$(conda shell.bash hook)"
conda activate "$VENV_DIR"
case "$stage" in
 prepare) python hierarchical_vamb_clustering_rocm.py --vamb-dir "$SAMPLE_DIR/vamb" --outdir "$OUTDIR" prepare-global ;;
 expand) python hierarchical_vamb_clustering_rocm.py --vamb-dir "$SAMPLE_DIR/vamb" --outdir "$OUTDIR" expand-global --shards 256 ;;
 fasta) python ~/atavide_lite/bin/vamb_create_fasta.py "vamb_crass/contigs.fna.gz" "$OUTDIR/hierarchical_unsplit.tsv" 20000 "$OUTDIR/hierarchical_unsplit_bins" ;;
 *) echo "Unknown stage: $stage" >&2; exit 2 ;;
esac
