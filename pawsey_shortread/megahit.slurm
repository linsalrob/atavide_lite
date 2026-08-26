#!/bin/bash
#SBATCH --job-name=MegaHit
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=230G


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

if [[ -s "megahit/$O/final.contigs.fa" && \
      -s "megahit/$O/final.graph.fastg" && \
      -s "megahit/$O/final.gfa" ]]; then
    echo "MEGAHIT outputs already exist for $O. Skipping." >&2
    exit 0
fi



#if [[ -e megahit/$O ]]; then
#	echo "megahit/$O already exists. Nothing to do" >&2;
#	exit 0;
#fi

echo "Using $HOSTREMOVED/$R1 $HOSTREMOVED/$R2 for the assembly" >&2;

#megahit -1 $HOSTREMOVED/$R1 -2 $HOSTREMOVED/$R2 -o megahit/$O -t $SLURM_CPUS_PER_TASK

#megahit_toolkit contig2fastg 119 megahit/$O/final.contigs.fa > megahit/$O/final.graph.fastg; 

#$HOME/GitHubs/atavide_lite/bin/fastg2gfa megahit/$O/final.graph.fastg > megahit/$O/final.gfa


# Run MEGAHIT only if final.contigs.fa is missing or empty.
if [[ -s "megahit/$O/final.contigs.fa" ]]; then
    echo "MEGAHIT contigs already exist for $O. Skipping assembly." >&2
elif [[ -d "megahit/$O" ]]; then
    echo "MEGAHIT directory exists for $O but final.contigs.fa is incomplete." >&2
    echo "Attempting to continue the existing assembly." >&2

    megahit \
        --continue \
        -o "megahit/$O"
else
    megahit \
        -1 "$HOSTREMOVED/$R1" \
        -2 "$HOSTREMOVED/$R2" \
        -o "megahit/$O" \
        -t "$SLURM_CPUS_PER_TASK"
fi

if [[ ! -s "megahit/$O/final.contigs.fa" ]]; then
    echo "ERROR: Missing or empty megahit/$O/final.contigs.fa" >&2
    exit 1
fi

# Generate FASTG only if it is missing or empty.
if [[ -s "megahit/$O/final.graph.fastg" ]]; then
    echo "FASTG already exists for $O. Skipping." >&2
else
    megahit_toolkit contig2fastg 119 \
        "megahit/$O/final.contigs.fa" \
        > "megahit/$O/final.graph.fastg"
fi

if [[ ! -s "megahit/$O/final.graph.fastg" ]]; then
    echo "ERROR: Missing or empty megahit/$O/final.graph.fastg" >&2
    exit 1
fi

# Generate GFA only if it is missing or empty.
if [[ -s "megahit/$O/final.gfa" ]]; then
    echo "GFA already exists for $O. Skipping." >&2
else
    "$HOME/GitHubs/atavide_lite/bin/fastg2gfa" \
        "megahit/$O/final.graph.fastg" \
        > "megahit/$O/final.gfa"
fi

if [[ ! -s "megahit/$O/final.gfa" ]]; then
    echo "ERROR: Missing or empty megahit/$O/final.gfa" >&2
    exit 1
fi

echo "MEGAHIT outputs completed successfully for $O." >&2
