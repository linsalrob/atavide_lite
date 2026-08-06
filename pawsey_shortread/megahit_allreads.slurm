#!/bin/bash
#SBATCH --job-name=MegaHitAll
#SBATCH --partition=highmem
#SBATCH --time=4-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=960G
#SBATCH -o slurm_output/%x-%j.out
#SBATCH -e slurm_output/%x-%j.err



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

find $HOSTREMOVED -name \*${FILEEND/R1/R2}

# Collect R1 and R2 lists
R1_FILES=($(find $HOSTREMOVED -name \*${FILEEND} | sort))
R2_FILES=($(find $HOSTREMOVED -name \*${FILEEND/R1/R2} | sort))
mapfile -t OTHER_FILES < <(
    comm -23 \
        <(find "$HOSTREMOVED" -name "*.fastq.gz" | sort) \
        <(find "$HOSTREMOVED" \( -name "*${FILEEND}" -o -name "*${FILEEND/R1/R2}" \) | sort)
)

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



