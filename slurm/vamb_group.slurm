#!/bin/bash
#SBATCH --job-name=VAMBgroups
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH -o vamb_groups-%j.out
#SBATCH -e vamb_groups-%j.err
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu

set -uo pipefail
eval "$(conda shell.bash hook)"
conda activate atavide_lite_vamb

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

	mkdir --parents $BGFS/$OUTDIR
	rsync $OUTDIR/contigs.fna.gz $BGFS/$OUTDIR/contigs.fna.gz 
	vamb -o C -p 16 --cuda --outdir $BGFS/$OUTDIR/vamb --fasta $BGFS/$OUTDIR/contigs.fna.gz --bamfiles $OUTDIR/mapped_reads/*bam 
	if [[ $? != 0 ]]; then 
		echo "Running VAMB on $OUTDIR failed." >&2;
	else
		python ~/atavide_lite/bin/vamb_create_fasta.py $BGFS/$OUTDIR/contigs.fna.gz $BGFS/$OUTDIR/vamb/vae_clusters.tsv 20000 $BGFS/$OUTDIR/bins
		rm -f $BGFS/$OUTDIR/contigs.fna.gz
		rsync -a $BGFS/$OUTDIR/ $OUTDIR/
	fi
done
