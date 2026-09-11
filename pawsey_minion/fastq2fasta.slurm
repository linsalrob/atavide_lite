#!/bin/bash
#SBATCH --job-name=fastq2fasta
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1800M
#SBATCH -o slurm_output/faing-%j.out
#SBATCH -e slurm_output/faing-%j.err

# NOTE: Setonix enforces MaxMemPerCPU (work=1840M, highmem=7900M). A bare --mem
# larger than cpus-per-task * MaxMemPerCPU does not fail - Slurm silently raises the
# CPU count to cover it, so the job queues for far more CPUs than it uses. Request
# memory with --mem-per-cpu so the two can never disagree.

set -euo pipefail
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# test we have atavide_lite conda installed, and then activate it
eval "$(conda shell.bash hook)"
"$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh"
ATAVIDE_CONDA="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite"
conda activate $ATAVIDE_CONDA

# Match conversion workers to allocated CPUs; each extract_R1 is single-threaded.
THREADS=${SLURM_CPUS_PER_TASK:-8}

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh


# this requires ~/GitHubs/atavide_lite/bin/fastq2fasta.c be compiled with make, and it has built in reading gzip support

extract_R1 () {
	FQ=$1;
	FB=`basename $FQ .gz`;
	FA=${FB/fastq/fasta}
	"$HOME/GitHubs/atavide_lite/bin/fastq2fasta" -n 1 "$FQ" - | gzip -c > "fasta/$FA.gz";
}

export -f extract_R1
mkdir -p fasta

find -L $HOSTREMOVED -type f -name *fastq.gz | parallel -j "$THREADS" extract_R1
