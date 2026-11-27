#!/bin/bash
#SBATCH --job-name=VAMB
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --gres=gpu:4
#SBATCH --partition=gpu



set -euo pipefail
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

THREADS=${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-8}}
echo "Slurm CPUs on node: ${SLURM_CPUS_ON_NODE:-not set}" >&2
echo "VAMB threads: $THREADS" >&2

# test we have atavide_lite conda installed, and then activate it
eval "$(conda shell.bash hook)"
"$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh"
ATAVIDE_CONDA="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite"
conda activate $ATAVIDE_CONDA

# Please see the vamb_install.slurm script to create the 
# virtual environment for vamb. NOTE: Do not use conda to install vamb, because 
# we need a rcom install, not nvidia.

# vamb_install.slurm builds this environment; keep the two paths in sync.
VENV_DIR="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/vamb"
if [[ ! -x "$VENV_DIR/bin/vamb" ]]; then
	echo "Please run vamb_install.slurm to create the vamb virtual environment" >&2
	exit 1;
fi

# The ROCm torch wheels bundle their own ROCm libraries, so no rocm module is needed
# (Setonix's 2026.11 stack dropped rocm/6.2.4). --stack keeps atavide_lite's minimap2 on PATH.
conda activate --stack "$VENV_DIR"

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh


OUTDIR=vamb
mkdir --parents $OUTDIR

#if [[ ! -e $OUTDIR/contigs.fna.gz ]]; then
#	python ~/GitHubs/atavide_lite/bin/vamb_concatenate.py -m 300 $OUTDIR/contigs.fna.gz megahit/*/final.contigs.fa 
#fi

if [[ ! -s "$OUTDIR/contigs.fna.gz" ]]; then
    echo "ERROR: $OUTDIR/contigs.fna.gz not found." >&2
    echo "Run vamb_concat.slurm before vamb.slurm." >&2
    exit 2
fi

mkdir --parents $OUTDIR $OUTDIR/mapped_reads  

if [ ! -e $OUTDIR/contigs.mmi ]; then minimap2 -I100G -d $OUTDIR/contigs.mmi $OUTDIR/contigs.fna.gz; fi

for R1 in $(sort -R R1_reads.txt); do
	R2=${R1/_R1/_R2}
	BAM=${R1/$FILEEND/.bam}
	echo "R1: $R2 BAM: $BAM" >&2;
	if [[ ! -e $OUTDIR/mapped_reads/$BAM ]]; then
		echo "Bam file: '$OUTDIR/mapped_reads/$BAM' not found. Please check if you want that file" >&2;
	fi
done



# VAMB refuses to write into an existing --outdir. If a previous run finished, keep
# its clusters; if it left a partial directory behind, stop rather than guess.
if [[ -s "$OUTDIR/vamb/vae_clusters_unsplit.tsv" && -s "$OUTDIR/vamb/vae_clusters_split.tsv" ]]; then
	echo "$OUTDIR/vamb has both cluster files. Skipping VAMB binning" >&2
elif [[ -e "$OUTDIR/vamb" ]]; then
	echo "ERROR: $OUTDIR/vamb exists but is missing vae_clusters_unsplit.tsv or vae_clusters_split.tsv (a failed run?)." >&2
	echo "       Remove it before re-running." >&2
	exit 2
else
	# Reuse a composition.npz saved from an earlier VAMB run (copy it to $OUTDIR/) to skip
	# recomputing tetranucleotide frequencies; otherwise compute them from the contigs.
	if [[ -s "$OUTDIR/composition.npz" ]]; then
		COMPOSITION=(--composition "$OUTDIR/composition.npz")
	else
		COMPOSITION=(--fasta "$OUTDIR/contigs.fna.gz")
	fi
	vamb bin default -o C -p "$THREADS" --cuda --outdir "$OUTDIR/vamb" "${COMPOSITION[@]}" --bamdir "$OUTDIR/mapped_reads/"
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

