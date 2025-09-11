#!/bin/bash
#SBATCH --job-name=MegaHitAll
#SBATCH --partition=highmem
#SBATCH --time=4-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=960G
#SBATCH -o slurm_output/megahit_slurm/%x-%j.out
#SBATCH -e slurm_output/megahit_slurm/%x-%j.err


set -euo pipefail
# set -eu
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

if [[ ! -n "$ATAVIDE_CONDA" ]]; then
    echo "Please define the ATAVIDE_CONDA environment variable with the location of your ATAVIDE_LITE conda installation" >&2;
    exit 2;
fi


if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh

eval "$(conda shell.bash hook)"
conda activate $ATAVIDE_CONDA




# Collect R1 and R2 lists
R1_FILES=($(find $HOSTREMOVED -name \*_R1.fastq.gz | sort))
R2_FILES=($(find $HOSTREMOVED -name \*_R2.fastq.gz | sort))
OTHER_FILES=($(find $HOSTREMOVED -name \*.fastq.gz | grep -v -E "_R[12]\.fastq\.gz" | sort)) || true

# Join them into comma-separated strings
R1_JOINED=$(IFS=,; echo "${R1_FILES[*]}")
R2_JOINED=$(IFS=,; echo "${R2_FILES[*]}")
OTHER_JOINED=$(IFS=,; echo "${OTHER_FILES[*]}")

echo -e "R1 files: ${R1_JOINED}\nR2 files: ${R2_JOINED}\nOthers: ${OTHER_JOINED}" >&2


# Build the final string

FINAL_STRING=""
if [[ ! -z  "${R1_JOINED}" ]]; then
    FINAL_STRING+="-1 ${R1_JOINED} "
fi
if [[ ! -z  "${R2_JOINED}" ]]; then
    FINAL_STRING+="-2 ${R2_JOINED} "
fi
if [[ ! -z  "${OTHER_JOINED}" ]]; then
    FINAL_STRING+="-r ${OTHER_JOINED} "
fi

echo "Running megahit with the following parameters: megahit $FINAL_STRING -o megahit_output" >&2

# megahit $FINAL_STRING --continue -o megahit_output 
megahit $FINAL_STRING -o megahit_all_reads



