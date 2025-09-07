#!/bin/bash
#SBATCH --job-name=VAMB
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH -o slurm_output/vamb_slurm/map_reads-%A_%a.out
#SBATCH -e slurm_output/vamb_slurm/map_reads-%A_%a.err

# Run the minimap section for each of the groups

set -euo pipefail

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
	echo "Running minimap on the contigs in $OUTDIR" >&2;

	R1=$(head -n $SLURM_ARRAY_TASK_ID R1_reads.txt | tail -n 1)
	R2=${R1/_R1/_R2}
	BAM=${R1/$FILEEND/.bam}

	if [[ ! -e $OUTDIR/mapped_reads/$BAM ]]; then
		mkdir --parents $OUTDIR $OUTDIR/mapped_reads 
		mkdir --parents $OUTDIR/mapped_reads

		if [[ ! -e $OUTDIR/contigs.mmi ]]; then
			minimap2 -I100G -d $OUTDIR/contigs.mmi $OUTDIR/contigs.fna.gz
		fi

		minimap2 -t 32 -N 5 -ax sr $OUTDIR/contigs.mmi --split-prefix mmsplit$$ $HOSTREMOVED/$R1 $HOSTREMOVED/$R2 \
			| samtools view -F 3584 -b --threads 16  | samtools sort -@ 16 -o $OUTDIR/mapped_reads/$BAM -
		samtools index $OUTDIR/mapped_reads/$BAM
	fi
done
