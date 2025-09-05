#!/bin/bash
#SBATCH --job-name=VAMBgroups
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH -o slurm_output/vamb_slurm/vamb_groups-%j.out
#SBATCH -e slurm_output/vamb_slurm/vamb_groups-%j.err
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu

set -uo pipefail

if [[ ! -n "$ATAVIDE_CONDA_VAMB" ]]; then
    echo "Please define the ATAVIDE_CONDA_VAMB environment variable with the location of your ATAVIDE_LITE_VAMB conda installation" >&2;
    echo "Note: this is based on the atavide_vamb.yaml and is different from the regular atavide conda" >&2;
    exit 2;
fi

eval "$(conda shell.bash hook)"
conda activate $ATAVIDE_CONDA_VAMB

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

	mkdir --parents $OUTDIR
	vamb -o C -p 16 --cuda --outdir $OUTDIR/vamb --fasta $OUTDIR/contigs.fna.gz --bamfiles $OUTDIR/mapped_reads/*bam 
	if [[ $? != 0 ]]; then 
		echo "Running VAMB on $OUTDIR failed." >&2;
	else
		python ~/atavide_lite/bin/vamb_create_fasta.py $OUTDIR/contigs.fna.gz $OUTDIR/vamb/vae_clusters.tsv 20000 $OUTDIR/bins
	fi
done
