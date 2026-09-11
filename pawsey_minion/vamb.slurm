#!/bin/bash
#SBATCH --job-name=VAMB
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH -o slurm_output/vamb_slurm/vamb-%j.out
#SBATCH -e slurm_output/vamb_slurm/vamb-%j.err
#SBATCH --account=pawsey1018-gpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
# Setonix allocates one 8-core chiplet per GPU, so 8 CPUs matches --gres=gpu:1.
# Without this SLURM_CPUS_PER_TASK is unset and THREADS silently falls back to 16.
#SBATCH --cpus-per-task=8
#SBATCH --time=60:00

set -euo pipefail
trap 's=$?; echo "$0: Error on line "$LINENO": $BASH_COMMAND"; exit $s' ERR

# test we have atavide_lite conda installed, and then activate it
eval "$(conda shell.bash hook)"
"$HOME/GitHubs/atavide_lite/pawsey_lib/check_atavide_lite_env.sh"
ATAVIDE_CONDA="/scratch/$PAWSEY_PROJECT/$USER/software/miniconda3/atavide_lite_vamb"
conda activate $ATAVIDE_CONDA

# Derive thread count from the Slurm allocation rather than hardcoding it.
THREADS=${SLURM_CPUS_PER_TASK:-8}

# module load pytorch/2.2.0-rocm5.7.3
module load pytorch/2.7.1-rocm6.3.3

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

mkdir -p $OUTDIR/mapped_reads

for READS in $(sort -R reads.txt); do
	BAM=${READS/$FILEEND/.bam}
	echo "READS: $READS BAM: $BAM" >&2;
	if [[ ! -e $OUTDIR/mapped_reads/$BAM ]]; then
		echo "Bam file: '$OUTDIR/mapped_reads/$BAM' not found. Data generated" >&2;
		mkdir -p $OUTDIR/mapped_reads
		minimap2 -t "$THREADS" -N 5 -ax map-ont $OUTDIR/contigs.mmi --split-prefix mmsplit$$ fastq/$READS | samtools view -F 3584 -b --threads "$THREADS"  | samtools sort -@ "$THREADS" -o $OUTDIR/mapped_reads/$BAM -
	fi
done


vamb -o C -p "$THREADS" --cuda --outdir $OUTDIR/vamb --fasta $OUTDIR/contigs.fna.gz --bamfiles $OUTDIR/mapped_reads/*bam

python "$HOME/GitHubs/atavide_lite/bin/vamb_create_fasta.py" $OUTDIR/contigs.fna.gz $OUTDIR/vamb/vae_clusters.tsv 20000 $OUTDIR/bins
