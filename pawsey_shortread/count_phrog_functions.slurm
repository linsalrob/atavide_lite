#!/bin/bash
#SBATCH --job-name=CountPHROG
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -o slurm_output/count_phrog-%j.out
#SBATCH -e slurm_output/count_phrog-%j.err

set -euo pipefail
trap 's=$?; echo "$0: Error on line $LINENO: $BASH_COMMAND" >&2; exit $s' ERR

if [[ ! -f DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SAMPLENAME" >&2
	exit 2
fi

source DEFINITIONS.sh
: "${SAMPLENAME:?SAMPLENAME must be defined in DEFINITIONS.sh}"
: "${PAWSEY_PROJECT:?PAWSEY_PROJECT must be set}"

PHROG_DIR=/scratch/$PAWSEY_PROJECT/$USER/Databases/phrog
PHROG_ANNOT="$PHROG_DIR/phrog_annot_v4.tsv.gz"
if [[ ! -s "$PHROG_ANNOT" ]]; then
	echo "PHROG annotation database is missing or empty: $PHROG_ANNOT" >&2
	echo "Run mmseqs_add_phrog_function.slurm to download and validate the PHROG database." >&2
	exit 2
fi
if ! gzip -t "$PHROG_ANNOT"; then
	echo "PHROG annotation database is not a valid gzip file: $PHROG_ANNOT" >&2
	echo "Remove or replace the damaged file, then rerun mmseqs_add_phrog_function.slurm." >&2
	exit 2
fi

OUTPUT_DIR=phrog_functions
if [[ -e "$OUTPUT_DIR" ]]; then
	echo "$OUTPUT_DIR already exists; refusing to overwrite or mix results" >&2
	exit 2
fi
TMP_OUTPUT="${OUTPUT_DIR}.tmp.${SLURM_JOB_ID:-$$}"
if [[ -e "$TMP_OUTPUT" ]]; then
	echo "Temporary output path already exists: $TMP_OUTPUT" >&2
	exit 2
fi
mkdir -p "$TMP_OUTPUT"
trap 'status=$?; if [[ $status -eq 0 ]]; then rm -rf -- "$TMP_OUTPUT"; else echo "Preserving temporary output at $TMP_OUTPUT (exit $status)" >&2; fi' EXIT

export PYTHONPATH="${PYTHONPATH:+${PYTHONPATH}:}$HOME/GitHubs/atavide_lite"
python "$HOME/GitHubs/atavide_lite/bin/count_phrog_functions.py" \
	-d mmseqs -o "$TMP_OUTPUT" -n "$SAMPLENAME" -a "$PHROG_ANNOT"
find "$TMP_OUTPUT" -type f -name '*.tsv' -exec pigz -- {} +
mv "$TMP_OUTPUT" "$OUTPUT_DIR"
trap - EXIT
