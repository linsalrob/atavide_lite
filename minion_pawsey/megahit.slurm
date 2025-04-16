#!/bin/bash
#SBATCH --job-name=MegaHit
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -o slurm_output/megahit_slurm/megahit-%A_%a.out
#SBATCH -e slurm_output/megahit_slurm/megahit-%A_%a.err


set -euo pipefail

if [[ ! -n "$ATAVIDE_CONDA" ]]; then
    echo "Please include the ATAVIDE_CONDA location on your slurm submission: --export=ATAVIDE_CONDA=$ATAVIDE_CONDA" >&2
    exit 2;
fi
eval "$(conda shell.bash hook)"
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


megahit -r $HOSTREMOVED/$READS -o megahit/$OUTPUT -t 8

megahit_toolkit contig2fastg 119 megahit/$OUTPUT/final.contigs.fa > megahit/$OUTPUT/final.graph.fastg; 
~/atavide_lite/bin/fastg2gfa megahit/$OUTPUT/final.graph.fastg > megahit/$OUTPUT/final.gfa


