#!/bin/bash


if [[ $# != 1 ]]; then
	echo "$0 <directory of no host files>" >&2;
	exit 1;
fi

SOURCE=$1
DEST=mmseqs
SCRIPT=~/GitHubs/atavide_lite/nci_pbs/fasta/mmseqs_easy_taxonomy.pbs

FAFILEEND="_R1.fasta.gz"

for R1 in $(find $SOURCE -type f -name \*_R1\* -printf "%f\n"); do
	OUT=${R1/$FAFILEEND/}
	if [[ -e mmseqs/$OUT/ ]]; then echo "$R1 Nothing to do!" >&2; exit 0; fi

	R2=${R1/_R1/_R2}

	if [[ -e $SOURCE/$R2 ]]; then
		echo "qsub -v R1=$SOURCE/$R1,R2=$SOURCE/$R2 $SCRIPT"
		qsub -v R1=$SOURCE/$R1,R2=$SOURCE/$R2 $SCRIPT
	else
		echo "qsub -v R1=$SOURCE/$R1 $SCRIPT"
		qsub -v R1=$SOURCE/$R1 $SCRIPT
	fi
done
