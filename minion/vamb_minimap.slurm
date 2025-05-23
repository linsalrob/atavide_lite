#!/bin/bash
#SBATCH --job-name=VAMB
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH -o slurm_output/vamb_slurm/map_reads-%A_%a.out
#SBATCH -e slurm_output/vamb_slurm/map_reads-%A_%a.err

# Just the minimap section of VAMB which does not require GPUs and we can do in parallel

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate atavide_lite

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh

OUTDIR=vamb

mkdir --parents $OUTDIR/mapped_reads

if [[ ! -e $OUTDIR/contigs.fna.gz ]]; then
	echo "Please create the contigs with vamb_concatenate" >&2;
	exit 2;
fi

if [[ ! -e $OUTDIR/contigs.mmi ]]; then
	echo "Please index the contigs with something like:" >&2;
	echo "\tminimap2 -I100G -d $OUTDIR/contigs.mmi $OUTDIR/contigs.fna.gz" >&2;
	exit 3;
fi


READS=$(head -n $SLURM_ARRAY_TASK_ID reads.txt | tail -n 1)
BAM=${READS/$FILEEND/.bam}

rsync -avz $OUTDIR/ $BGFS/$OUTDIR/

if [[ ! -e $OUTDIR/mapped_reads/$BAM ]]; then
	mkdir --parents $BGFS/$OUTDIR/mapped_reads
	if [[ ! -e $BGFS/$OUTDIR/contigs.mmi ]]; then
		echo "Crap. $BGFS/$OUTDIR/contigs.mmi not found";
		exit 4;
	fi

	echo "minimap2 -t 16 -N 5 -ax map-ont $OUTDIR/contigs.mmi --split-prefix mmsplit$$ fastq/$READS | samtools view -F 3584 -b --threads 16  | samtools sort -@ 16 -o $OUTDIR/mapped_reads/$BAM -" >&2
	minimap2 -t 16 -N 5 -ax map-ont $BGFS/$OUTDIR/contigs.mmi --split-prefix $BGFS/mmsplit$$ fastq/$READS | samtools view -F 3584 -b --threads 16  | samtools sort -@ 16 -o $BGFS/$OUTDIR/mapped_reads/$BAM -
	samtools index $BGFS/$OUTDIR/mapped_reads/$BAM
fi

rsync -avz $BGFS/$OUTDIR/mapped_reads/ $OUTDIR/mapped_reads/
