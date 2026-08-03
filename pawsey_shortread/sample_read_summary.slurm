#!/bin/bash
#SBATCH --account=pawsey1018
#SBATCH --job-name=sample_summary
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32GB
#SBATCH -o slurm_output/sample_summary-%j.out
#SBATCH -e slurm_output/sample_summary-%j.err

set -euo pipefail
trap 's=$?; echo "$0: Error on line $LINENO: $BASH_COMMAND" >&2; exit $s' ERR

if [[ ! -f R1_reads.txt ]]; then
	echo "R1_reads.txt does not exist" >&2
	exit 2
fi

if [[ ! -f DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SAMPLENAME, FILEEND, and HOSTREMOVED" >&2
	exit 2
fi

source DEFINITIONS.sh
: "${SAMPLENAME:?SAMPLENAME must be defined in DEFINITIONS.sh}"

export PYTHONPATH="${PYTHONPATH:-}:$HOME/GitHubs/atavide_lite"

METADATA_ARGS=()
if [[ -f sample_metadata.tsv ]]; then
	METADATA_ARGS=(--metadata sample_metadata.tsv)
fi

python "$HOME/GitHubs/atavide_lite/bin/sample_read_summary.py" \
	--reads R1_reads.txt \
	--definitions DEFINITIONS.sh \
	--threads "$SLURM_CPUS_PER_TASK" \
	--output "${SAMPLENAME}_sample_read_summary.tsv" \
	"${METADATA_ARGS[@]}" \
	--verbose
