#!/bin/bash


if [[ $# != 1 ]]; then
	echo "$0 <directory of fasta files>" >&2;
	exit 1;
fi

SOURCE=$1
DEST=human
SCRIPT=~/GitHubs/atavide_lite/pbs/fasta/human.pbs
for FILE in $(find $SOURCE -name \*_R1* -printf "%f\n"); do 
	if [[ ! -e $DEST/$FILE ]]; then
		qsub -v R1=$SOURCE/$FILE $SCRIPT
	fi
done
