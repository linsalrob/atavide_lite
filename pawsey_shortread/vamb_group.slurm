#!/bin/bash
#SBATCH --job-name=VAMBgroups
#SBATCH --time=1-0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -o slurm_output/vamb_slurm/vamb_groups-%j.out
#SBATCH -e slurm_output/vamb_slurm/vamb_groups-%j.err
#SBATCH --gres=gpu:8
#SBATCH --partition=gpu

set -uo pipefail

# Please see the vamb_install.slurm script to create the 
# virtual environment for vamb. NOTE: Do not use conda, because 
# we need a rcom install, not nvidia.

VENV_DIR="/scratch/pawsey1018/edwa0468/software/pip/vamb"
if [[ ! -e "$VENV_DIR" ]]; then 
	echo "Please run vamb_install.slurm to create the vamb virtual environment" >&2
	exit 1;
fi

module load rocm/6.2.4
source "$VENV_DIR/bin/activate"

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh


VAMBDIR=vamb_groups

for INF in $(find $VAMBDIR -name contigs.fna.gz); do
	OUTDIR=`dirname $INF`;
	echo "Running vamb on $OUTDIR" >&2;

	if [[ ! -e $OUTDIR/contigs.fna.gz ]]; then
		echo "No contigs.fna.gz file found for $OUTDIR, which is unlikely!" >&2;
		continue;
	fi


	if [[ ! -e $OUTDIR/mapped_reads ]]; then
		echo "No mapped_reads found for $OUTDIR. Did you run vamb_minimap_group.slurm" >&2;
		continue;
	fi


	vamb bin default -o C -p 16 --cuda --outdir $OUTDIR/vamb --fasta $OUTDIR/contigs.fna.gz --bamdir $OUTDIR/mapped_reads/
	if [[ $? != 0 ]]; then 
		echo "Running VAMB on $OUTDIR failed." >&2;
	else
		python ~/atavide_lite/bin/vamb_create_fasta.py $OUTDIR/contigs.fna.gz $OUTDIR/vamb/vae_clusters.tsv 1 $OUTDIR/bins
	fi
done
