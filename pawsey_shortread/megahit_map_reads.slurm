#!/bin/bash
#SBATCH --job-name=MegaHitMap
#SBATCH --time=0-8
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH -o slurm_output/megahit_slurm/megahit-map-%A_%a.out
#SBATCH -e slurm_output/megahit_slurm/megahit-map-%A_%a.err

set -euo pipefail
trap 's=$?; echo "$0: Error on line $LINENO: $BASH_COMMAND" >&2; exit $s' ERR

# Test we have atavide_lite conda installed, and then activate it.
eval "$(conda shell.bash hook)"
"$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh"
ATAVIDE_CONDA="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite"
conda activate "$ATAVIDE_CONDA"

if [[ ! -f R1_reads.txt ]]; then
	echo "Please make a file with the R1 reads using this command:" >&2
	echo 'find fastq -name \*R1\* -printf "%f\n" > R1_reads.txt' >&2
	exit 2
fi

if [[ ! -f DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with FILEEND and HOSTREMOVED" >&2
	exit 2
fi

source DEFINITIONS.sh
: "${FILEEND:?FILEEND must be defined in DEFINITIONS.sh}"
: "${HOSTREMOVED:?HOSTREMOVED must be defined in DEFINITIONS.sh}"
: "${SLURM_ARRAY_TASK_ID:?Submit this script as a SLURM array job}"
: "${SLURM_CPUS_PER_TASK:?SLURM_CPUS_PER_TASK is not set}"

R1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" R1_reads.txt)
if [[ -z "$R1" ]]; then
	echo "No sample was found on line $SLURM_ARRAY_TASK_ID of R1_reads.txt" >&2
	exit 2
fi
R2=${R1/_R1/_R2}
SAMPLE=${R1%"$FILEEND"}
if [[ "$SAMPLE" == "$R1" ]]; then
	echo "ERROR: $R1 does not end with FILEEND ($FILEEND)" >&2
	exit 2
fi

ASSEMBLY_DIR="megahit/$SAMPLE"
REFERENCE="$ASSEMBLY_DIR/final.contigs.fa"
BAM="$ASSEMBLY_DIR/mapped_reads.bam"
TMP_BAM="$ASSEMBLY_DIR/.mapped_reads.${SLURM_JOB_ID:-$$}.bam"
TMP_PREFIX="${TMPDIR:-/tmp}/megahit-map-${SLURM_JOB_ID:-$$}-${SLURM_ARRAY_TASK_ID}"

for INPUT in "$HOSTREMOVED/$R1" "$HOSTREMOVED/$R2" "$REFERENCE"; do
	if [[ ! -s "$INPUT" ]]; then
		echo "Missing or empty input: $INPUT" >&2
		exit 2
	fi
done

if [[ -s "$BAM" && -s "$BAM.bai" ]] && samtools quickcheck "$BAM"; then
	echo "$BAM and its index already exist and passed samtools quickcheck. Nothing to do." >&2
	exit 0
fi
if [[ -e "$BAM" || -e "$BAM.bai" ]]; then
	echo "Incomplete output exists at $BAM or $BAM.bai; refusing to overwrite it" >&2
	exit 2
fi

cleanup() {
	status=$?
	rm -f -- "$TMP_BAM" "$TMP_BAM.bai"
	exit "$status"
}
trap cleanup EXIT

echo "Mapping $HOSTREMOVED/$R1 and $HOSTREMOVED/$R2 to $REFERENCE" >&2
minimap2 -t "$SLURM_CPUS_PER_TASK" --secondary=no -ax sr \
	--split-prefix "$TMP_PREFIX" "$REFERENCE" "$HOSTREMOVED/$R1" "$HOSTREMOVED/$R2" \
	| samtools view -@ "$SLURM_CPUS_PER_TASK" -b -F 4 - \
	| samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$TMP_BAM" -

samtools index -@ "$SLURM_CPUS_PER_TASK" "$TMP_BAM"
samtools quickcheck "$TMP_BAM"
mv -- "$TMP_BAM" "$BAM"
mv -- "$TMP_BAM.bai" "$BAM.bai"

echo "Created $BAM and $BAM.bai" >&2
