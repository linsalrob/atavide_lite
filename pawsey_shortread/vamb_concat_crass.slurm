#!/bin/bash
#SBATCH --job-name=VambConcat
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH -o slurm_output/vamb_slurm/%x-%j.out
#SBATCH -e slurm_output/vamb_slurm/%x-%j.err

# we use cross assembly of all samples using megahit, and then map the reads on a per group
# basis

set -euo pipefail

VENV_DIR="/scratch/$PAWSEY_ACCOUNT/$USER/software/pip/vamb"
if [[ ! -e "$VENV_DIR" ]]; then 
	echo "Please run vamb_install.slurm to create the vamb virtual environment" >&2
	exit 1;
fi

module load rocm/6.2.4
source "$VENV_DIR/bin/activate"


if [[ ! -e R1_reads.txt ]]; then 
	echo "Please make a file with the R1 reads using this command:" >&2
	echo "find fastq -name \*R1\* -printf "%f\n" > R1_reads.txt" >&2;
	exit 2;
fi


if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh

eval "$(conda shell.bash hook)"
conda activate $ATAVIDE_CONDA

# run vamb_concat 
# we want a tsv file with [category] [sample]
# category should be one word, no spaces, and should provide a single grep. We will test this
# and not continue.
# sample should be unique in the sample

E_BADARGS=65

MEGAHITDIR="megahit_all_reads"
VAMB="vamb_crass"

if [[ "$#" -ne 1 ]]; then
	echo "$0 <file with [group, sample]>";
	exit $E_BADARGS;
fi

if [[ ! -d $MEGAHITDIR ]]; then
	echo "Can't find $MEGAHITDIR. Please run megahit_host_removed.sh first" >&2;
	exit 2;
fi

if [[ ! -e $MEGAHITDIR/final.contigs.fa ]]; then
	echo "Can't find $MEGAHITDIR/final.contigs.fa. Did megahit finish successfully?" >&2;
	exit 2;
fi


IFS=$'\n'

declare -A FQPATHS

mkdir -p $VAMB

for LN in $(grep -v Sample $1); do
	GRP=$(echo $LN | cut -d$'\t' -f 1)
	SMP=$(echo $LN | cut -d$'\t' -f 2)
	# check how many reads we have for $SMP. Should be only 1
	LC=0;
	LC=$(grep $SMP R1_reads.txt | wc -l) || echo "No matches";
	if [[ $LC == 0 ]]; then
		echo "ERROR: Can't find $SMP in R1_reads.txt. Skipped" >&2;
		continue;
	fi
	if [[ $LC -gt 1 ]]; then
		echo "ERROR: More than one R1 file for $SMP. (We found $LC files in R1_reads.txt)" >&2;
		echo "We can't continue to process this data. Can you please check $SMP and see if there ">&2;
		echo " is a way to make it unique in this file (e.g. add _ or -)" >&2;
		exit 2;
	fi

	mkdir -p $VAMB/$GRP
	echo "python ~/atavide_lite/bin/vamb_concatenate.py -m 300 $VAMB/$GRP/contigs.fna.gz $MEGAHITDIR/final.contigs.fa" >&2;
	bash -c "python ~/atavide_lite/bin/vamb_concatenate.py -m 300 $VAMB/$GRP/contigs.fna.gz $MEGAHITDIR/final.contigs.fa";
done


