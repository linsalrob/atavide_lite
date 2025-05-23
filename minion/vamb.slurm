#!/bin/bash
#SBATCH --job-name=VAMB
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH -o slurm_output/vamb_slurm/vamb-%j.out
#SBATCH -e slurm_output/vamb_slurm/vamb-%j.err
#SBATCH --gres=gpu:1
#SBATCH --time=60:00

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate atavide_lite_vamb


if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh


OUTDIR=vamb
mkdir --parents $OUTDIR

if [[ ! -e $OUTDIR/contigs.fna.gz ]]; then
	echo "Please run vamb concatenate" 2>&1;
	exit 1;
fi

for READS in $(sort -R reads.txt); do
	BAM=${READS/$FILEEND/.bam}
	echo "READS: $READS BAM: $BAM" >&2;
	if [[ -e $OUTDIR/mapped_reads/$BAM ]]; then
		echo "Bam file: '$OUTDIR/mapped_reads/$BAM' not found. Data generated" >&2;
		minimap2 -t 16 -N 5 -ax map-ont $OUTDIR/contigs.mmi --split-prefix mmsplit$$ fastq/$READS | samtools view -F 3584 -b --threads 16  | samtools sort -@ 16 -o $OUTDIR/mapped_reads/$BAM -
	fi
done


cp $OUTDIR/contigs.fna.gz $BGFS/
cp $OUTDIR/mapped_reads/*bam $BGFS/

vamb -o C -p 16 --cuda --outdir $BGFS/$OUTDIR/vamb --fasta $BGFS/contigs.fna.gz --bamfiles $BGFS/*bam 

python ~/atavide_lite/bin/vamb_create_fasta.py $BGFS/contigs.fna.gz $BGFS/$OUTDIR/vamb/vae_clusters.tsv 20000 $BGFS/$OUTDIR/bins

rsynv -avz $BGFS/$OUTDIR/ $OUTDIR/
