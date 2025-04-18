#!/bin/bash
#SBATCH --job-name=CountSS
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=short
#SBATCH -o count_ss-%j.out
#SBATCH -e count_ss-%j.err


if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh

python ~/atavide_lite/bin/count_subsystems.py -d mmseqs -s subsystems -n $SAMPLENAME
find subsystems -type f -exec pigz {} \;
