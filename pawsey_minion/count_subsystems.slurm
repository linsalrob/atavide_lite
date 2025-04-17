#!/bin/bash
#SBATCH --job-name=CountSS
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH -o slurm_output/count_ss-%j.out
#SBATCH -e slurm_output/count_ss-%j.err


if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh

perl ~/atavide_lite/bin/count_subsystems.pl -d mmseqs -n $SAMPLENAME
