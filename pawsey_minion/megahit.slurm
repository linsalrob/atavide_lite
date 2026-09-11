#!/bin/bash
#SBATCH --job-name=MegaHit
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem-per-cpu=1800M
#SBATCH -o slurm_output/megahit_slurm/megahit-%A_%a.out
#SBATCH -e slurm_output/megahit_slurm/megahit-%A_%a.err

# NOTE: Setonix enforces MaxMemPerCPU (work=1840M). A bare --mem larger than
# cpus-per-task * MaxMemPerCPU does not fail - Slurm silently raises the CPU count to
# cover it, which also changes $SLURM_CPUS_PER_TASK and therefore THREADS below.

set -euo pipefail
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# test we have atavide_lite conda installed, and then activate it
eval "$(conda shell.bash hook)"
"$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh"
ATAVIDE_CONDA="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite"
conda activate $ATAVIDE_CONDA

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh

mkdir --parents megahit
READS=$(head -n $SLURM_ARRAY_TASK_ID reads.txt | tail -n 1)
OUTPUT=${READS/$FILEEND/};


if [[ -e megahit/$OUTPUT ]]; then
	echo "megahit/$OUTPUT already exists. Nothing to do" >&2;
	exit 0;
fi

echo "Using $HOSTREMOVED/$READS for the assembly" >&2;


THREADS=${SLURM_CPUS_PER_TASK:-8}
megahit -r $HOSTREMOVED/$READS -o megahit/$OUTPUT -t "$THREADS"

megahit_toolkit contig2fastg 119 megahit/$OUTPUT/final.contigs.fa > megahit/$OUTPUT/final.graph.fastg; 
"$HOME/GitHubs/atavide_lite/bin/fastg2gfa" megahit/$OUTPUT/final.graph.fastg > megahit/$OUTPUT/final.gfa

