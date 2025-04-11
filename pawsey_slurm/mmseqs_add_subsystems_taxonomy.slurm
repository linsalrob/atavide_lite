#!/bin/bash
#SBATCH --job-name=mmseqs_ss
#SBATCH --time=0-10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -o slurm_output/mmseqs_slurm/mmseqs_ss-%A_%a.out
#SBATCH -e slurm_output/mmseqs_slurm/mmseqs_ss-%A_%a.err


if [[ ! -n "$ATAVIDE_CONDA" ]]; then
    echo "Please define the ATAVIDE_CONDA environment variable with the location of your ATAVIDE_LITE conda installation" >&2;
    exit 2;
fi

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

# if we are job #1 we are going to copy the database, otherwise we
# are going to wait for the "copying_done" flag

DBDIR=~/scratch_1018/uniref
eval "$(conda shell.bash hook)"

if [[ $SLURM_ARRAY_TASK_ID == 1 ]]; then
	if [ ! -e $DBDIR/copying_done ]; then
		conda activate rclone
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

conda activate $ATAVIDE_CONDA


R1=$(head -n $SLURM_ARRAY_TASK_ID R1_reads.txt | tail -n 1)
BASE=${R1/$FILEEND/}

echo "Running ~/atavide_lite/bin/easy_taxonomy_to_function_taxa.py -f mmseqs/$BASE/${BASE}_tophit_report.gz -d $BGFS/uniref.sqlite > $BGFS/${BASE}_tophit_report_subsystems_taxa" >&2
python ~/atavide_lite/bin/easy_taxonomy_to_function_taxa.py -f mmseqs/$BASE/${BASE}_tophit_report.gz -d $BGFS/uniref.sqlite > $BGFS/${BASE}_tophit_report_subsystems_taxa
pigz $BGFS/${BASE}_tophit_report_subsystems_taxa 
mkdir -p mmseqs/$BASE
mv -f $BGFS/${BASE}_tophit_report_subsystems_taxa mmseqs/$BASE/${BASE}_tophit_report_subsystems_taxa.gz

