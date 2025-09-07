#!/bin/bash
#SBATCH --job-name=VAMB
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH -o slurm_output/vamb_slurm/vamb-%j.out
#SBATCH -e slurm_output/vamb_slurm/vamb-%j.err
#SBATCH --gres=gpu:8
#SBATCH --partition=gpu

set -euo pipefail

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


OUTDIR=vamb
mkdir --parents $OUTDIR

if [[ ! -e $OUTDIR/contigs.fna.gz ]]; then
	python ~/atavide_lite/bin/vamb_concatenate.py -m 300 $OUTDIR/contigs.fna.gz megahit/*/output/final.contigs.fa 
fi

mkdir --parents $OUTDIR $OUTDIR/mapped_reads  

if [ ! -e $OUTDIR/contigs.mmi ]; then minimap2 -I100G -d $OUTDIR/contigs.mmi $OUTDIR/contigs.fna.gz; fi

for R1 in $(sort -R R1_reads.txt); do
	R2=${R1/_R1/_R2}
	BAM=${R1/$FILEEND/.bam}
	echo "R1: $R2 BAM: $BAM" >&2;
	if [[ ! -e $OUTDIR/mapped_reads/$BAM ]]; then
		echo "Bam file: '$OUTDIR/mapped_reads/$BAM' not found. Data generated" >&2;
		minimap2 -t 16 -N 5 -ax sr $OUTDIR/contigs.mmi --split-prefix mmsplit$$ $HOSTREMOVED/$R1 $HOSTREMOVED/$R2 | samtools view -F 3584 -b --threads 16  | samtools sort -@ 16 -o $OUTDIR/mapped_reads/$BAM -
	fi
done



vamb -o C -p 16 --cuda --outdir $OUTDIR/vamb --fasta $OUTDIR/contigs.fna.gz --bamfiles $OUTDIR/mapped_reads/*bam 

python ~/atavide_lite/bin/vamb_create_fasta.py $OUTDIR/contigs.fna.gz $OUTDIR/vamb/vae_clusters.tsv 20000 $OUTDIR/bins

