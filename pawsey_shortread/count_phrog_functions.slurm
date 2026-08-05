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
TMP_OUTPUT="${OUTPUT_DIR}.tmp.${SLURM_JOB_ID:-$$}"
if [[ -e "$TMP_OUTPUT" ]]; then
	echo "Temporary output path already exists: $TMP_OUTPUT" >&2
	exit 2
fi
mkdir -- "$TMP_OUTPUT"
PUBLISH_LOCK="${OUTPUT_DIR}.publish.lock"
PUBLISH_LOCK_HELD=0
cleanup() {
	status=$?
	if [[ $PUBLISH_LOCK_HELD -eq 1 ]]; then
		if ! rmdir -- "$PUBLISH_LOCK"; then
			echo "Unable to remove publish lock $PUBLISH_LOCK" >&2
			[[ $status -ne 0 ]] || status=1
		fi
	fi
	if [[ $status -eq 0 ]]; then
		rm -rf -- "$TMP_OUTPUT"
	else
		echo "Preserving temporary output at $TMP_OUTPUT (exit $status)" >&2
	fi
	exit "$status"
}
trap cleanup EXIT

export PYTHONPATH="${PYTHONPATH:+${PYTHONPATH}:}$HOME/GitHubs/atavide_lite"
python "$HOME/GitHubs/atavide_lite/bin/count_phrog_functions.py" \
	-d mmseqs -o "$TMP_OUTPUT" -n "$SAMPLENAME" -a "$PHROG_ANNOT"
find "$TMP_OUTPUT" -type f -name '*.tsv' -exec pigz -- {} +

# mkdir is the exclusive check-and-publish lock: only one concurrent job can
# inspect OUTPUT_DIR and rename its completed temporary directory at a time.
if ! mkdir -- "$PUBLISH_LOCK"; then
	echo "Another job is publishing $OUTPUT_DIR; preserving $TMP_OUTPUT" >&2
	exit 2
fi
PUBLISH_LOCK_HELD=1
if [[ -e "$OUTPUT_DIR" ]]; then
	echo "$OUTPUT_DIR already exists; refusing to overwrite or mix results" >&2
	exit 2
fi
mv "$TMP_OUTPUT" "$OUTPUT_DIR"
