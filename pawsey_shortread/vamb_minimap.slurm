#!/bin/bash
#SBATCH --job-name=V_minimap
#SBATCH --time=0-8
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G


# Just the minimap section of VAMB which does not require GPUs and we can do in parallel


set -euo pipefail
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# test we have atavide_lite conda installed, and then activate it
eval "$(conda shell.bash hook)"
"$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh"
ATAVIDE_CONDA="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite"
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


R1=$(head -n $SLURM_ARRAY_TASK_ID R1_reads.txt | tail -n 1)
R2=${R1/_R1/_R2}
BAM=${R1/$FILEEND/.bam}

if [[ ! -s "$OUTDIR/mapped_reads/$BAM" ]]; then
    mkdir --parents "$OUTDIR/mapped_reads"

    if [[ ! -s "$OUTDIR/contigs.mmi" ]]; then
        echo "ERROR: Minimap2 index is missing: $OUTDIR/contigs.mmi" >&2
        exit 2
    fi

    minimap2 \
        -t "$SLURM_CPUS_PER_TASK" \
        --secondary=no \
        --split-prefix "$OUTDIR/mmsplit${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}" \
        -ax sr \
        "$OUTDIR/contigs.mmi" \
        "$HOSTREMOVED/$R1" \
        "$HOSTREMOVED/$R2" |
        samtools view \
            -b \
            -F 3596 \
            --threads "$SLURM_CPUS_PER_TASK" |
        samtools sort \
            -@ "$SLURM_CPUS_PER_TASK" \
            -o "$OUTDIR/mapped_reads/$BAM" \
            -

    samtools index "$OUTDIR/mapped_reads/$BAM"
fi
