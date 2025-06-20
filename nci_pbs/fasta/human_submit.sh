#!/bin/bash


if [[ $# != 1 ]]; then
	echo "$0 <directory of fasta files>" >&2;
	exit 1;
fi

SOURCE=$1
DEST=human
SCRIPT=~/GitHubs/atavide_lite/nci_pbs/fasta/human.pbs
for FILE in $(find $SOURCE -type f -name \*_R1\* -printf "%f\n"); do 
	if [[ -e $SOURCE/${FILE/R1/R2} ]] && [[ ! -e $DEST/$FILE ]]; then
		qsub -v R1=$SOURCE/$FILE,R2=$SOURCE/${FILE/R1/R2} $SCRIPT
	elif [[ ! -e $DEST/$FILE ]]; then
		qsub -v R1=$SOURCE/$FILE $SCRIPT
	fi
done
