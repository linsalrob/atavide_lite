#!/bin/bash
#SBATCH --job-name=mmseqs_phrog
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -o slurm_output/mmseqs_slurm/%x-%A_%a.out
#SBATCH -e slurm_output/mmseqs_slurm/%x-%A_%a.err

set -euo pipefail
trap 's=$?; echo "$0: Error on line $LINENO: $BASH_COMMAND" >&2; exit $s' ERR

PHROG_DIR=/scratch/$PAWSEY_PROJECT/$USER/Databases/phrog
PHROG_ANNOT="$PHROG_DIR/phrog_annot_v4.tsv.gz"
PHROG_UNIPROT="$PHROG_DIR/phrog_uniprot.tsv.gz"

# Download either database file if it is absent. The lock prevents array tasks
# from trying to populate the shared directory at the same time.
if [[ ! -f "$PHROG_ANNOT" || ! -f "$PHROG_UNIPROT" ]]; then
	eval "$(conda shell.bash hook)"
	mkdir -p "$PHROG_DIR"

	exec 9>"$PHROG_DIR/.phrog_download.lock"
	flock 9
	if [[ ! -f "$PHROG_ANNOT" ]]; then
		conda run -n rclone -- rclone copy A18:dbs/phrog/phrog_annot_v4.tsv.gz "$PHROG_DIR/"
	fi
	if [[ ! -f "$PHROG_UNIPROT" ]]; then
		conda run -n rclone -- rclone copy A18:dbs/phrog/phrog_uniprot.tsv.gz "$PHROG_DIR/"
	fi
	flock -u 9
fi

for database in "$PHROG_ANNOT" "$PHROG_UNIPROT"; do
	if [[ ! -s "$database" ]] || ! gzip -t "$database"; then
		echo "$database is missing, empty, or not a valid gzip file" >&2
		exit 1
	fi
done

if [[ ! -f R1_reads.txt ]]; then
	echo "R1_reads.txt does not exist" >&2
	exit 2
fi

if [[ ! -f DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file defining FILEEND" >&2
	exit 2
fi

source DEFINITIONS.sh
: "${FILEEND:?FILEEND must be defined in DEFINITIONS.sh}"
: "${SLURM_ARRAY_TASK_ID:?Submit this script as a SLURM array job}"

R1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" R1_reads.txt)
if [[ -z "$R1" ]]; then
	echo "No sample was found on line $SLURM_ARRAY_TASK_ID of R1_reads.txt" >&2
	exit 2
fi

BASE=${R1%"$FILEEND"}
INPUT="mmseqs/$BASE/${BASE}_tophit_report.gz"
OUTPUT="mmseqs/$BASE/${BASE}_tophit_report_phrog_function.tsv.gz"

if [[ ! -f "$INPUT" ]]; then
	echo "$INPUT does not exist" >&2
	exit 2
fi

if [[ -e "$OUTPUT" ]]; then
	echo "$OUTPUT already exists. Not overwriting" >&2
	exit 0
fi

TMP_OUTPUT="${OUTPUT}.tmp.${SLURM_JOB_ID:-$$}.${SLURM_ARRAY_TASK_ID}"
trap 'rm -f "$TMP_OUTPUT"' EXIT

echo "$(date): Adding PHROG annotations to $INPUT" >&2
{
	gzip -cd "$PHROG_ANNOT"
	printf '__PHROG_SECTION__\n'
	gzip -cd "$PHROG_UNIPROT"
	printf '__MMSEQS_SECTION__\n'
	gzip -cd "$INPUT"
} | awk -F '\t' -v OFS='\t' '
	$0 == "__PHROG_SECTION__" { section = 1; next }
	$0 == "__MMSEQS_SECTION__" { section = 2; next }
	section == 0 {
		if ($1 != "phrog") annotation[$1] = $2 OFS $3 OFS $4
		next
	}
	section == 1 {
		phrog_number = $1
		sub(/^phrog_/, "", phrog_number)
		# Keep the first entry if a UniRef ID occurs more than once.
		if (!($2 in uniref_phrog)) {
			uniref_phrog[$2] = $1
			uniref_annotation[$2] = annotation[phrog_number]
		}
		next
	}
	{
		if ($1 in uniref_phrog) {
			print $0, uniref_phrog[$1], uniref_annotation[$1]
		} else {
			print $0, "", "", "", ""
		}
	}
' | gzip -c > "$TMP_OUTPUT"

gzip -t "$TMP_OUTPUT"
mv "$TMP_OUTPUT" "$OUTPUT"
trap - EXIT
echo "$(date): Wrote $OUTPUT" >&2
