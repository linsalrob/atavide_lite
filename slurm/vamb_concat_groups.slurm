#!/bin/bash
#SBATCH --job-name=VambConcat
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=short
#SBATCH -o vamb_concat-%j.out
#SBATCH -e vamb_concat-%j.err

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate atavide_lite_vamb

# run vamb_concat 
# we want a tsv file with [category] [sample]
# category should be one word, no spaces, and should provide a single grep. We will test this
# and not continue.
# sample should be unique in the sample

E_BADARGS=65

if [[ "$#" -ne 1 ]]; then
	echo "$0 <file with [group, sample]>";
	exit $E_BADARGS;
fi

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED" >&2
	exit 2;
fi

source DEFINITIONS.sh


IFS=$'\n'

declare -A FQPATHS

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
	
	R1=$(grep $SMP R1_reads.txt)
	FQ=${R1/$FILEEND/};

	
	if [[ ! -e "megahit/$FQ/output/final.contigs.fa" ]]; then
		echo "No contigs file, megahit/$FQ/output/final.contigs.fa" >&2;
		continue
	fi

	if [[ "${FQPATHS[$GRP]:-UNBOUND}" != "UNBOUND" ]]; then
		FQPATHS[$GRP]+=" megahit/$FQ/output/final.contigs.fa";
	else
		#${FQPATHS[$GRP]} = "megahit/$FQ/output/final.contigs.fa";
		FQPATHS+=([$GRP]="megahit/$FQ/output/final.contigs.fa");
	fi
done


mkdir -p vamb_groups
for LC in ${!FQPATHS[@]}; do 
	echo "Making $LC" >&2;
	mkdir -p vamb_groups/$LC
	echo "python ~/atavide_lite/bin/vamb_concatenate.py vamb_groups/$LC/contigs.fna.gz ${FQPATHS[$LC]}" >&2;
	bash -c "python ~/atavide_lite/bin/vamb_concatenate.py vamb_groups/$LC/contigs.fna.gz ${FQPATHS[$LC]}";
done


