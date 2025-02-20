#!/bin/bash
#SBATCH --job-name=mmseqs_ss
#SBATCH --time=7-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -o slurm_output/mmseqs_slurm/mmseqs_ss-%A_%a.out
#SBATCH -e slurm_output/mmseqs_slurm/mmseqs_ss-%A_%a.err

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


R1=$(head -n $SLURM_ARRAY_TASK_ID R1_reads.txt | tail -n 1)
BASE=${R1/$FILEEND/}

if [[ -e mmseqs/$BASE/${BASE}_tophit_report_subsystems_taxa.gz ]]; then 
	echo "mmseqs/$BASE/${BASE}_tophit_report_subsystems_taxa.gz already exists. Not overwriting" >&2;
	exit 2
fi

cp /home/edwa0468/UniRef/uniref.sqlite $BGFS
cp mmseqs/$BASE/${BASE}_tophit_report.gz $BGFS

echo `date` "Running ~/atavide_lite/bin/easy_taxonomy_to_function_taxa.py -f mmseqs/$BASE/${BASE}_tophit_report.gz -d $BGFS/uniref.sqlite > $BGFS/${BASE}_tophit_report_subsystems_taxa" >&2
python ~/atavide_lite/bin/easy_taxonomy_to_function_taxa.py -f $BGFS/${BASE}_tophit_report.gz -d $BGFS/uniref.sqlite > $BGFS/${BASE}_tophit_report_subsystems_taxa
echo `date` "easy_taxonomy_to_function_taxa complete" >&2

pigz $BGFS/${BASE}_tophit_report_subsystems_taxa 
mkdir -p mmseqs/$BASE
mv -f $BGFS/${BASE}_tophit_report_subsystems_taxa.gz mmseqs/$BASE/${BASE}_tophit_report_subsystems_taxa.gz

