#!/bin/bash
#SBATCH --job-name=mmseqs_ssf
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH -o slurm_output/mmseqs_slurm/mmseqs_ss-%A_%a.out
#SBATCH -e slurm_output/mmseqs_slurm/mmseqs_ss-%A_%a.err

set -euo pipefail


eval "$(conda shell.bash hook)"
conda activate atavide_lite


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

DBDIR=$HOME/UniRef/

cp $DBDIR/uniref.sqlite $BGFS/

R1=$(head -n $SLURM_ARRAY_TASK_ID reads.txt | tail -n 1)
BASE=${R1/$FILEEND/}

if [[ -e mmseqs/$BASE/${BASE}_tophit_report_subsystems_taxa.gz ]]; then 
	echo "mmseqs/$BASE/${BASE}_tophit_report_subsystems_taxa.gz already exists. Not overwriting" >&2;
	exit 2
fi

mkdir -p $BGFS/mmseqs/$BASE/
cp mmseqs/$BASE/${BASE}_tophit_report.gz $BGFS/mmseqs/$BASE/
echo `date` "Running ~/atavide_lite/bin/easy_taxonomy_to_function_taxa_fast.py -f $BGFS/mmseqs/$BASE/${BASE}_tophit_report.gz -d $BGFS/uniref.sqlite -o $BGFS/mmseqs/$BASE/${BASE}_tophit_report_subsystems_taxa" >&2
python ~/atavide_lite/bin/easy_taxonomy_to_function_taxa_fast.py -f $BGFS/mmseqs/$BASE/${BASE}_tophit_report.gz -d $BGFS/uniref.sqlite -o $BGFS/mmseqs/$BASE/${BASE}_tophit_report_subsystems_taxa
echo `date` "easy_taxonomy_to_function_taxa complete" >&2

pigz $BGFS/mmseqs/$BASE/${BASE}_tophit_report_subsystems_taxa 
mv $BGFS/mmseqs/$BASE/${BASE}_tophit_report_subsystems_taxa.gz mmseqs/$BASE/${BASE}_tophit_report_subsystems_taxa
