#!/bin/bash
#SBATCH --job-name=MegaHit
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH -o slurm_output/megahit_slurm/megahit-%A_%a.out
#SBATCH -e slurm_output/megahit_slurm/megahit-%A_%a.err

set -euo pipefail
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# test we have atavide_lite conda installed, and then activate it
eval "$(conda shell.bash hook)"
"$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh"
ATAVIDE_CONDA="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite"
conda activate $ATAVIDE_CONDA

if [[ ! -e R1_reads.txt ]]; then 
	echo "Please make a file with the R1 reads using this command:" >&2
	echo "find fastq -name \*R1\* -printf "%f\n" > R1_reads.txt" >&2;
	exit 2;
fi


if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh

mkdir --parents megahit
R1=$(head -n $SLURM_ARRAY_TASK_ID R1_reads.txt | tail -n 1)
R2=${R1/_R1/_R2}
O=${R1/$FILEEND/};


if [[ -e megahit/$O ]]; then
	echo "megahit/$O already exists. Nothing to do" >&2;
	exit 0;
fi

echo "Using $HOSTREMOVED/$R1 $HOSTREMOVED/$R2 for the assembly" >&2;

megahit -1 $HOSTREMOVED/$R1 -2 $HOSTREMOVED/$R2 -o megahit/$O -t 36

megahit_toolkit contig2fastg 119 megahit/$O/final.contigs.fa > megahit/$O/final.graph.fastg; 
~/atavide_lite/bin/fastg2gfa megahit/$O/final.graph.fastg > megahit/$O/final.gfa


