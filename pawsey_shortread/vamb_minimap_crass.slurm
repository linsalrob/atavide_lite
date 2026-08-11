#!/bin/bash
#SBATCH --job-name=VAMB
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH -o slurm_output/vamb_slurm/map_reads-%A_%a.out
#SBATCH -e slurm_output/vamb_slurm/map_reads-%A_%a.err

# Run the minimap section for each of the groups

# Submit this as:
# sbatch --array=1-$NUM_R1_READS:1 --account=${PAWSEY_PROJECT} $PAWSEY_SRC/vamb_minimap_crass.slurm samples.tsv
# Only runs minimap2 so does not need GPU nodes.


set -euo pipefail
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR


SAMPLES=$1
VAMBDIR=vamb_crass

if [[ ! -e $SAMPLES ]]; then
	echo "Please provide a samples.tsv file as the first argument" >&2
	echo "We expect a tab-delimited file with [group]\t[R1 reads] format" >&2
	echo "The group should be a non-whitespace string as we make a directory in the $VAMBDIR for each group" >&2
	echo "The R1 reads should be the name of the R1 fastq file (e.g. sample_R1.fastq.gz)" >&2
	exit 2;
fi

if [[ ! -e $VAMBDIR ]]; then
	echo "$VAMBDIR does not exist. Please check which directory you want to process" >&2;
	exit 2;
fi

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


# Before we begin, lets index the contigs.fna.gz file with minimap2. We will do this only once, and then all the other jobs will wait for it to finish.
INDEX="$VAMBDIR/contigs.mmi"
INDEX_DONE="$VAMBDIR/contigs.mmi.done"
INDEX_FAILED="$VAMBDIR/contigs.mmi.failed"

if [[ "$SLURM_ARRAY_TASK_ID" == "$SLURM_ARRAY_TASK_MIN" ]]; then

    rm -f "$INDEX_FAILED"

    if [[ ! -e "$INDEX" ]]; then
        echo "[$(date)] Building minimap2 index"

        if ! minimap2 -I100G \
             -d "${INDEX}.tmp" \
             "$VAMBDIR/contigs.fna.gz"
        then
            rm -f "${INDEX}.tmp"
            touch "$INDEX_FAILED"
            exit 1
        fi

        mv "${INDEX}.tmp" "$INDEX"
    fi

    touch "$INDEX_DONE"

else

    echo "[$(date)] Waiting for index"

    while [[ ! -e "$INDEX_DONE" ]]; do
        if [[ -e "$INDEX_FAILED" ]]; then
            echo "ERROR: minimap2 indexing failed"
            exit 1
        fi
        sleep 10
    done

fi

# Here, we are going to read each line of the samples.tsv file.
# They have [group]\t[R1 reads] format.
# If the group directory does not exist, we create it and then
# use minimap to make the bam files in that directory.
#


LN=$(head -n $SLURM_ARRAY_TASK_ID $SAMPLES | tail -n 1)


GRP=$(echo $LN | cut -d$'\t' -f 1)
SMP=$(echo $LN | cut -d$'\t' -f 2)
OUTDIR=$VAMBDIR/$GRP
R2=${R1/_R1/_R2}
BAM=${R1/$FILEEND/.bam}

if [[ ! -e $OUTDIR/mapped_reads/$BAM ]]; then
	mkdir --parents $OUTDIR/mapped_reads

	minimap2 -t 32 -N 5 -ax sr $VAMBDIR/contigs.mmi --split-prefix mmsplit$$ $HOSTREMOVED/$R1 $HOSTREMOVED/$R2 \
		| samtools view -b -F 3584 --threads 16  | samtools sort -@ 16 -o $OUTDIR/mapped_reads/$BAM -
	samtools index $OUTDIR/mapped_reads/$BAM
else
	echo "Skipped $OUTDIR/mapped_reads/$BAM as it already exists" >&2;
fi
