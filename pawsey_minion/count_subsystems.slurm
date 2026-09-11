#!/bin/bash
#SBATCH --job-name=CountSS
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1800M
#SBATCH -o slurm_output/count_ss-%j.out
#SBATCH -e slurm_output/count_ss-%j.err

# NOTE: Setonix enforces MaxMemPerCPU (work=1840M, highmem=7900M). A bare --mem
# larger than cpus-per-task * MaxMemPerCPU does not fail - Slurm silently raises the
# CPU count to cover it, so the job queues for far more CPUs than it uses. Request
# memory with --mem-per-cpu so the two can never disagree.


if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh

export PYTHONPATH=$PYTHONPATH:$HOME/GitHubs/atavide_lite
python "$HOME/GitHubs/atavide_lite/bin/count_subsystems.py" -d mmseqs -n $SAMPLENAME
find subsystems -type f -exec pigz {} \;
