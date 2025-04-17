#!/bin/bash
#SBATCH --job-name=mmseqs_ss
#SBATCH --time=0-10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -o slurm_output/mmseqs_slurm/mmseqs_ss-%A_%a.out
#SBATCH -e slurm_output/mmseqs_slurm/mmseqs_ss-%A_%a.err

if [[ ! -n "$ATAVIDE_CONDA" ]]; then
    echo "Please include the ATAVIDE_CONDA location on your slurm submission: --export=ATAVIDE_CONDA=$ATAVIDE_CONDA" >&2
    exit 2;
fi
eval "$(conda shell.bash hook)"
conda activate $ATAVIDE_CONDA


if [[ ! -e reads.txt ]]; then 
	echo "Please make a file with the reads using this command:" >&2
	echo "find fastq -type f -printf "%f\n" > reads.txt" >&2;
	exit 2;
fi


if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh

# if we are job #1 we are going to copy the database, otherwise we
# are going to wait for the "copying_done" flag

DBDIR=~/scratch_1018/uniref

if [[ $SLURM_ARRAY_TASK_ID == 1 ]]; then
	if [ ! -e $DBDIR/copying_done ]; then
		mkdir -p $DBDIR
		rclone copy A18:dbs/uniref.sqlite $DBDIR/
		touch $DBDIR/copying_done
	fi
else
	while [ ! -e $DBDIR/copying_done ]; do
		# sleep a random amount of time between 1 and 60 seconds. 
		# Hopefully this will also cause random offsets in reading the db
		sleep $(( (RANDOM % 30) + 1 ))
	done
fi


READS=$(head -n $SLURM_ARRAY_TASK_ID reads.txt | tail -n 1)
BASE=${READS/$FILEEND/}

echo "Running ~/atavide_lite/bin/easy_taxonomy_to_function_taxa.py -f mmseqs/$BASE/${BASE}_tophit_report.gz -d $DBDIR/uniref.sqlite > mmseqs/$BASE/${BASE}_tophit_report_subsystems_taxa" >&2
python ~/atavide_lite/bin/easy_taxonomy_to_function_taxa.py -f mmseqs/$BASE/${BASE}_tophit_report.gz -d $DBDIR/uniref.sqlite > mmseqs/$BASE/${BASE}_tophit_report_subsystems_taxa
pigz mmseqs/$BASE//${BASE}_tophit_report_subsystems_taxa 

