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

THREADS=${SLURM_CPUS_PER_TASK:-16}

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh


OUTDIR=vamb
mkdir --parents $OUTDIR/mapped_reads

if [[ ! -e $OUTDIR/contigs.fna.gz ]]; then
	echo "Please run vamb concatenate" >&2;
	exit 1;
fi

# Normally vamb_minimap.slurm has already made every BAM; this only fills gaps.
# Map the host-removed reads: they are what megahit assembled, and mapping raw reads
# would let host reads inflate contig coverage.
for READS in $(sort -R reads.txt); do
	BAM=${READS/$FILEEND/.bam}
	echo "READS: $READS BAM: $BAM" >&2;
	if [[ ! -e $OUTDIR/mapped_reads/$BAM ]]; then
		echo "Bam file: '$OUTDIR/mapped_reads/$BAM' not found. Data generated" >&2;
		if [[ ! -s $OUTDIR/contigs.mmi ]]; then
			minimap2 -d $OUTDIR/contigs.mmi $OUTDIR/contigs.fna.gz
		fi
		minimap2 -t $THREADS -N 5 -ax map-ont $OUTDIR/contigs.mmi --split-prefix $BGFS/mmsplit$$ "$HOSTREMOVED/$READS" | samtools view -F 3584 -b --threads $THREADS  | samtools sort -@ $THREADS -o $OUTDIR/mapped_reads/$BAM -
	fi
done

# VAMB >= 4 takes "vamb bin default"; the bare "vamb --fasta ..." form no longer parses.
# It refuses to write into an existing --outdir, so skip a finished run (both cluster files
# present) and stop on a partial one. Only the VAMB step is skipped: the bins below are
# always rebuilt from the cluster files, so a run interrupted while copying back cannot
# leave them missing.
if [[ -s "$OUTDIR/vamb/vae_clusters_unsplit.tsv" && -s "$OUTDIR/vamb/vae_clusters_split.tsv" ]]; then
	echo "$OUTDIR/vamb has both cluster files. Skipping VAMB binning" >&2
elif [[ -e "$OUTDIR/vamb" ]]; then
	echo "ERROR: $OUTDIR/vamb exists but is missing vae_clusters_unsplit.tsv or vae_clusters_split.tsv" >&2
	echo "       (a failed or interrupted run?). Remove it before re-running." >&2
	exit 2
else
	mkdir --parents $BGFS/$OUTDIR/mapped_reads
	cp $OUTDIR/contigs.fna.gz $BGFS/$OUTDIR/
	cp $OUTDIR/mapped_reads/*.bam $BGFS/$OUTDIR/mapped_reads/

	vamb bin default -o C -p $THREADS --cuda --outdir $BGFS/$OUTDIR/vamb --fasta $BGFS/$OUTDIR/contigs.fna.gz --bamdir $BGFS/$OUTDIR/mapped_reads/

	# The BAMs are already in $OUTDIR/mapped_reads; copy back only what VAMB produced.
	rsync -a $BGFS/$OUTDIR/vamb $OUTDIR/
fi

# Split bins (one per sample) go to bins/, which checkm.slurm reads; unsplit bins are kept
# apart so CheckM does not score every contig twice.
python ~/GitHubs/atavide_lite/bin/vamb_create_fasta.py $OUTDIR/contigs.fna.gz $OUTDIR/vamb/vae_clusters_unsplit.tsv 20000 $OUTDIR/bins_unsplit &
UNSPLIT_PID=$!
python ~/GitHubs/atavide_lite/bin/vamb_create_fasta.py $OUTDIR/contigs.fna.gz $OUTDIR/vamb/vae_clusters_split.tsv 20000 $OUTDIR/bins &
SPLIT_PID=$!
# a bare "wait" returns 0 even when a job failed, so wait on each one
wait "$UNSPLIT_PID"
wait "$SPLIT_PID"
