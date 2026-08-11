#!/bin/bash
#SBATCH --job-name=VAMBcrass
#SBATCH --time=1-0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -o slurm_output/vamb_slurm/%x-%A_%a.out
#SBATCH -e slurm_output/vamb_slurm/%x-%$_%a.err
#SBATCH --gres=gpu:8
#SBATCH --partition=gpu

set -uo pipefail

# Please see the vamb_install.slurm script to create the 
# virtual environment for vamb. NOTE: Do not use conda, because 
# we need a rcom install, not nvidia.
#
# Please submit this script as an array job with the number of samples in the vamb_crass directory.
# sbatch --array=1-$(find vamb_crass -maxdepth 1 -mindepth 1 -type d | wc -l) --account=${PAWSEY_PROJECT}-gpu $PAWSEY_SRC/vamb_crass.slurm

VENV_DIR="/scratch/$PAWSEY_PROJECT/$USER/software/pip/vamb"
if [[ ! -e "$VENV_DIR" ]]; then 
	echo "Please run vamb_install.slurm to create the vamb virtual environment" >&2
	exit 1;
fi

module load rocm/6.2.4
source "$VENV_DIR/bin/activate"

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh


VAMBDIR=vamb_crass

if [[ ! -e $VAMBDIR ]]; then
	echo "Please run vamb_minimap_crass.slurm to create the $VAMBDIR directory" >&2
	exit 3;
fi

if [[ ! -e $VAMBDIR/contigs.fna.gz ]]; then
	echo "Please run vamb_minimap_crass.slurm to create the $VAMBDIR/contigs.fna.gz file" >&2
	exit 4;
fi

# find the top level directories
INPUTS=($(find  $VAMBDIR -maxdepth 1 -mindepth 1 -type d))
SAMPLE=${INPUTS[$SLURM_ARRAY_TASK_ID-1]}

if [[ ! -e $SAMPLE/mapped_reads ]]; then
	echo "No mapped_reads found for $SAMPLE. Did you run vamb_minimap_crass.slurm" >&2;
	exit 5;
fi

echo "Running vamb on $SAMPLE" >&2;

# check for empty bamfiles and remove them from mapped_reads
for BAM in $(find --bamdir $SAMPLE/mapped_reads/ -name \*bam); do 
	C=$(samtools idxstats $BAM | awk '{s+=$3} END{print s}'); 
	if [[ $C == 0 ]]; then
		EMPDIR=$(dirname $BAM); 
		mkdir -p $SAMPLE/empty_bamfiles/$EMPDIR; 
		echo "Moving empty bamfile $BAM to $SAMPLE/empty_bamfiles/$EMPDIR" >&2;
		mv -t $SAMPLE/empty_bamfiles/$EMPDIR $BAM; 
	fi;
done


vamb bin default -o C -p 16 --cuda --outdir $SAMPLE/vamb --fasta $VAMBDIR/contigs.fna.gz --bamdir $SAMPLE/mapped_reads/ -m 250
if [[ $? != 0 ]]; then 
	echo "Running VAMB on $SAMPLE failed." >&2;
else
	python ~/atavide_lite/bin/vamb_create_fasta.py $VAMBDIR/contigs.fna.gz $SAMPLE/vamb/vae_clusters_unsplit.tsv 20000 $SAMPLE/unsplit_bins &
	python ~/atavide_lite/bin/vamb_create_fasta.py $VAMBDIR/contigs.fna.gz $SAMPLE/vamb/vae_clusters_split.tsv 20000 $SAMPLE/split_bins &
	wait
fi
