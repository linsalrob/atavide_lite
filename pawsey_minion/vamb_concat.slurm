#!/bin/bash
#SBATCH --job-name=VambConcat
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1800M
#SBATCH -o slurm_output/vamb_slurm/vamb_concat-%j.out
#SBATCH -e slurm_output/vamb_slurm/vamb_concat-%j.err

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

# vamb_concatenate.py imports vamb, which lives in the environment built by
# pawsey_shortread/vamb_install.slurm.
VENV_DIR="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/vamb"
if [[ ! -x "$VENV_DIR/bin/vamb" ]]; then
	echo "Please run pawsey_shortread/vamb_install.slurm to create the VAMB environment at $VENV_DIR" >&2
	exit 1
fi
conda activate --stack "$VENV_DIR"

mkdir -p vamb
python "$HOME/GitHubs/atavide_lite/bin/vamb_concatenate.py" -m 300 vamb/contigs.fna.gz megahit/*/final.contigs.fa
