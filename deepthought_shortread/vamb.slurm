#!/bin/bash
#SBATCH --job-name=VAMB
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH -o slurm_output/vamb_slurm/vamb-%j.out
#SBATCH -e slurm_output/vamb_slurm/vamb-%j.err
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate atavide_lite

THREADS=${SLURM_CPUS_PER_TASK:-16}

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh


OUTDIR=vamb
mkdir --parents $OUTDIR

if [[ ! -e $OUTDIR/contigs.fna.gz ]]; then
	python ~/GitHubs/atavide_lite/bin/vamb_concatenate.py -m 300 $OUTDIR/contigs.fna.gz megahit/*/output/final.contigs.fa
fi

mkdir --parents $BGFS/$OUTDIR $BGFS/$OUTDIR/mapped_reads  $BGFS/fastq
mkdir --parents $OUTDIR/mapped_reads

cp $OUTDIR/contigs.fna.gz $BGFS/$OUTDIR
minimap2 -I100G -d $BGFS/$OUTDIR/contigs.mmi $BGFS/$OUTDIR/contigs.fna.gz
for R1 in $(sort -R R1_reads.txt); do
	R2=${R1/_R1/_R2}
	BAM=${R1/$FILEEND/.bam}
	echo "R1: $R2 BAM: $BAM" >&2;
	if [[ -e $OUTDIR/mapped_reads/$BAM ]]; then
		echo "Bam file: $BAM found. Data copied" >&2;
		cp $OUTDIR/mapped_reads/$BAM $BGFS/$OUTDIR/mapped_reads/$BAM;
	else
		echo "Bam file: '$OUTDIR/mapped_reads/$BAM' not found. Data generated" >&2;
		cp $HOSTREMOVED/$R1 $HOSTREMOVED/$R2 $BGFS/fastq
		minimap2 -t $THREADS -N 5 -ax sr $BGFS/$OUTDIR/contigs.mmi --split-prefix $BGFS/mmsplit$$ $BGFS/fastq/$R1 $BGFS/fastq/$R2 | samtools view -F 3584 -b --threads $THREADS  | samtools sort -@ $THREADS -o $BGFS/$OUTDIR/mapped_reads/$BAM -
	fi
	rsync -a $BGFS/$OUTDIR/ $OUTDIR/
done

rsync -a $BGFS/$OUTDIR/ $OUTDIR/

# conda activate atavide_lite_vamb
conda activate vamb

# VAMB >= 4 takes "vamb bin default"; the bare "vamb --fasta ..." form no longer parses.
# It refuses to write into an existing --outdir, so keep a finished run and stop on a
# partial one.
if [[ -s "$OUTDIR/vamb/vae_clusters_unsplit.tsv" && -s "$OUTDIR/vamb/vae_clusters_split.tsv" ]]; then
	echo "$OUTDIR/vamb has both cluster files. Skipping VAMB binning" >&2
elif [[ -e "$OUTDIR/vamb" ]]; then
	echo "ERROR: $OUTDIR/vamb exists but is missing vae_clusters_unsplit.tsv or vae_clusters_split.tsv (a failed run?)." >&2
	echo "       Remove it before re-running." >&2
	exit 2
else
	vamb bin default -o C -p $THREADS --cuda --outdir $OUTDIR/vamb --fasta $BGFS/$OUTDIR/contigs.fna.gz --bamdir $BGFS/$OUTDIR/mapped_reads/
fi

# Split bins (one per sample) go to bins/, which checkm.slurm reads; unsplit bins are kept
# apart so CheckM does not score every contig twice.
python ~/GitHubs/atavide_lite/bin/vamb_create_fasta.py $BGFS/$OUTDIR/contigs.fna.gz $OUTDIR/vamb/vae_clusters_unsplit.tsv 20000 $OUTDIR/bins_unsplit &
UNSPLIT_PID=$!
python ~/GitHubs/atavide_lite/bin/vamb_create_fasta.py $BGFS/$OUTDIR/contigs.fna.gz $OUTDIR/vamb/vae_clusters_split.tsv 20000 $OUTDIR/bins &
SPLIT_PID=$!
# a bare "wait" returns 0 even when a job failed, so wait on each one
wait "$UNSPLIT_PID"
wait "$SPLIT_PID"

rsync -a $BGFS/$OUTDIR/ $OUTDIR/
