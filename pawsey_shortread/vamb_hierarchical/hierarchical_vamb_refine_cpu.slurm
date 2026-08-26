#!/bin/bash
# CPU preparation/merge/FASTA stages for CheckM-guided refinement.
#SBATCH --job-name=hierVAMBrefCPU
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
sample=${1:?Usage: hierarchical_vamb_refine_cpu.slurm SAMPLE [prepare|merge|fasta]}
stage=${2:?Specify prepare, merge, or fasta}
VENV_DIR="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/vamb"
SAMPLE_DIR="vamb_crass/$sample"
OUTDIR="$SAMPLE_DIR/hierarchical_vamb"
QA_PATH="${REFINEMENT_QA:-$OUTDIR/checkm/qa.tsv}"
REFINEMENT_DIR="${REFINEMENT_DIR:-marker_refinement}"
eval "$(conda shell.bash hook)"
conda activate "$VENV_DIR"
case "$stage" in
 prepare) python refine_hierarchical_vamb.py --vamb-dir "$SAMPLE_DIR/vamb" --outdir "$OUTDIR" --refinement-dir "$REFINEMENT_DIR" prepare --qa "$QA_PATH" --assignments "$OUTDIR/hierarchical_unsplit.tsv" ;;
 merge) python refine_hierarchical_vamb.py --vamb-dir "$SAMPLE_DIR/vamb" --outdir "$OUTDIR" --refinement-dir "$REFINEMENT_DIR" merge --assignments "$OUTDIR/hierarchical_unsplit.tsv" ;;
 fasta) python ~/atavide_lite/bin/vamb_create_fasta.py "vamb_crass/contigs.fna.gz" "$OUTDIR/$REFINEMENT_DIR/refined_unsplit.tsv" 20000 "$OUTDIR/$REFINEMENT_DIR/refined_unsplit_bins" ;;
 *) echo "Unknown stage: $stage" >&2; exit 2 ;;
esac
