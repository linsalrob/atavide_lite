#!/bin/bash
#SBATCH --job-name=VambConcat
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=230G


set -euo pipefail

# test we have atavide_lite conda installed, and then activate it
eval "$(conda shell.bash hook)"
"$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh"
ATAVIDE_CONDA="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite"
conda activate $ATAVIDE_CONDA

VENV_DIR="/scratch/$PAWSEY_PROJECT/$USER/software/pip/vamb"

if [[ ! -e "$VENV_DIR" ]]; then
    echo "Please run vamb_install.slurm to create the VAMB virtual environment" >&2
    exit 1
fi

module load rocm/6.2.4
source "$VENV_DIR/bin/activate"

OUTDIR="vamb"
mkdir -p "$OUTDIR"

if [[ -s "$OUTDIR/contigs.fna.gz" ]]; then
    echo "$OUTDIR/contigs.fna.gz already exists. Skipping VAMB concatenation." >&2
else
    python "$HOME/GitHubs/atavide_lite/bin/vamb_concatenate.py" \
        -m 2000 \
        "$OUTDIR/contigs.fna.gz" \
        megahit/*/final.contigs.fa
fi

if [[ ! -s "$OUTDIR/contigs.fna.gz" ]]; then
    echo "ERROR: Failed to generate $OUTDIR/contigs.fna.gz" >&2
    exit 2
fi

if [[ ! -s "$OUTDIR/contigs.mmi" ||
      "$OUTDIR/contigs.fna.gz" -nt "$OUTDIR/contigs.mmi" ]]; then

    echo "Building minimap2 index" >&2

    rm -f "$OUTDIR/contigs.mmi.tmp"

    minimap2 \
        -I12G \
        -d "$OUTDIR/contigs.mmi.tmp" \
        "$OUTDIR/contigs.fna.gz"

    mv "$OUTDIR/contigs.mmi.tmp" "$OUTDIR/contigs.mmi"
else
    echo "$OUTDIR/contigs.mmi already exists and is current." >&2
fi
