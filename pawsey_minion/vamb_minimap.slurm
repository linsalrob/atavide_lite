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
if [[ ! -n "$ATAVIDE_CONDA" ]]; then
    echo "Please include the ATAVIDE_CONDA location on your slurm submission: --export=ATAVIDE_CONDA=$ATAVIDE_CONDA" >&2
    exit 2;
fi
eval "$(conda shell.bash hook)"
conda activate $ATAVIDE_CONDA


if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh

OUTDIR=vamb
mkdir --parents $OUTDIR

if [[ ! -e $OUTDIR/contigs.fna.gz ]]; then
	echo "Please create the contigs with vamb_concatenate" >&2;
	exit 2;
fi


READS=$(head -n $SLURM_ARRAY_TASK_ID reads.txt | tail -n 1)
BAM=${READS/$FILEEND/.bam}

if [[ ! -e $OUTDIR/mapped_reads/$BAM ]]; then
	mkdir --parents $OUTDIR $OUTDIR/mapped_reads 
	mkdir --parents $OUTDIR/mapped_reads

	if [[ ! -e $OUTDIR/contigs.mmi ]]; then
		minimap2 -I100G -d $OUTDIR/contigs.mmi $OUTDIR/contigs.fna.gz
	fi

	minimap2 -t 16 -N 5 -ax map-ont $OUTDIR/contigs.mmi --split-prefix mmsplit$$ fastq/$READS | samtools view -F 3584 -b --threads 16  | samtools sort -@ 16 -o $OUTDIR/mapped_reads/$BAM -
	samtools index $OUTDIR/mapped_reads/$BAM
fi
