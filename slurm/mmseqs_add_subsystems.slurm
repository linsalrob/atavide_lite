#!/bin/bash
#SBATCH --job-name=mmseqs_ss
#SBATCH --time=0-10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH -o slurm_output/mmseqs_slurm/mmseqs_ss-%A_%a.out
#SBATCH -e slurm_output/mmseqs_slurm/mmseqs_ss-%A_%a.err
#SBATCH --partition=short

eval "$(conda shell.bash hook)"
conda activate atavide_lite

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


cp /home/edwa0468/UniRef/uniref.sqlite $BGFS

R1=$(head -n $SLURM_ARRAY_TASK_ID R1_reads.txt | tail -n 1)
BASE=${R1/$FILEEND/}

python ~/atavide_lite/bin/easy_taxonomy_to_function.py -f mmseqs/$BASE/${BASE}_tophit_report.gz -d $BGFS/uniref.sqlite | gzip -c >  mmseqs/$BASE/${BASE}_tophit_report_subsystems.gz
