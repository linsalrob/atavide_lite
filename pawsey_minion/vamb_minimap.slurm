#!/bin/bash
#SBATCH --job-name=VAMB
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1800M
#SBATCH -o slurm_output/vamb_slurm/map_reads-%A_%a.out
#SBATCH -e slurm_output/vamb_slurm/map_reads-%A_%a.err

# NOTE: Setonix enforces MaxMemPerCPU (work=1840M, highmem=7900M). A bare --mem
# larger than cpus-per-task * MaxMemPerCPU does not fail - Slurm silently raises the
# CPU count to cover it, so the job queues for far more CPUs than it uses. Request
# memory with --mem-per-cpu so the two can never disagree.

# Just the minimap section of VAMB which does not require GPUs and we can do in parallel

set -euo pipefail
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# test we have atavide_lite conda installed, and then activate it
eval "$(conda shell.bash hook)"
"$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh"
ATAVIDE_CONDA="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite"
conda activate $ATAVIDE_CONDA

# Derive thread count from the Slurm allocation rather than hardcoding it.
THREADS=${SLURM_CPUS_PER_TASK:-16}

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

	minimap2 -t "$THREADS" -N 5 -ax map-ont $OUTDIR/contigs.mmi --split-prefix mmsplit$$ fastq/$READS | samtools view -F 3584 -b --threads "$THREADS"  | samtools sort -@ "$THREADS" -o $OUTDIR/mapped_reads/$BAM -
	samtools index $OUTDIR/mapped_reads/$BAM
fi
