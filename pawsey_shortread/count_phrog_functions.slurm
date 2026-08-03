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

export PYTHONPATH="${PYTHONPATH:+${PYTHONPATH}:}$HOME/GitHubs/atavide_lite"
python "$HOME/atavide_lite/bin/count_phrog_functions.py" -d mmseqs -n "$SAMPLENAME"
find phrog_functions -type f -name '*.tsv' -exec pigz {} +
