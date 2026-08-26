#!/bin/bash
# CheckM wrapper for hierarchical VAMB bins. Large bin sets can create a
# placement tree deeper than Python's default recursion limit.

#SBATCH --job-name=checkm
#SBATCH --time=2-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB
#SBATCH --partition=highmem
#SBATCH -o slurm_output/checkm-%j.out
#SBATCH -e slurm_output/checkm-%j.err

set -euo pipefail
trap 's=$?; echo "$0: Error on line $LINENO: $BASH_COMMAND"; exit $s' ERR

eval "$(conda shell.bash hook)"
"$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh"
ATAVIDE_CONDA="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite"
conda activate "$ATAVIDE_CONDA"

if [[ -z ${CHECKM_DATA_PATH:-} ]]; then
    CHECKM_DATA_PATH="$ATAVIDE_CONDA/checkm_data"
    [[ -e "$CHECKM_DATA_PATH" ]] || {
        echo "Missing CheckM data: $CHECKM_DATA_PATH" >&2
        exit 1
    }
fi

if [[ $# -lt 2 ]]; then
    echo "$0 <source directory> <destination directory>" >&2
    exit 1
fi

SOURCE=${1%/}
REPORT_DIR=${2%/}
SITE_DIR="${SLURM_TMPDIR:-/tmp}/checkm-sitecustomize-${SLURM_JOB_ID}"
mkdir -p "$SITE_DIR"
printf '%s\n' 'import sys' 'sys.setrecursionlimit(max(sys.getrecursionlimit(), 50000))' \
    > "$SITE_DIR/sitecustomize.py"
export PYTHONPATH="$SITE_DIR${PYTHONPATH:+:$PYTHONPATH}"
WORK_DIR="${SLURM_TMPDIR:-/tmp}/checkm-${SLURM_JOB_ID}"
mkdir -p "$WORK_DIR/tmp" "$REPORT_DIR"
threads=16

echo "Running CheckM on $SOURCE; local work: $WORK_DIR; report: $REPORT_DIR; threads: $threads" >&2
# CheckM writes many per-bin HMM files. Keep them on node-local storage so
# Lustre's per-user inode quota is not consumed; retain only the final QA TSV.
checkm lineage_wf --tmpdir "$WORK_DIR/tmp" -t "$threads" "$SOURCE" "$WORK_DIR"
checkm qa --tab_table -o 1 -f "$REPORT_DIR/qa.tsv" --tmpdir "$WORK_DIR/tmp" \
    -t "$threads" "$WORK_DIR/lineage.ms" "$WORK_DIR"
