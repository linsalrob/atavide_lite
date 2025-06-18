#!/bin/bash


if [[ $# != 1 ]]; then
	echo "$0 <directory of (split) fasta files>" >&2;
	exit 1;
fi

SOURCE=$1
DEST=cutadapt
SCRIPT=~/GitHubs/atavide_lite/pbs/fasta/cutadapt.pbs

for FILE in $(find $SOURCE -type f -printf "%f\n"); do 
	if [[ ! -e $DEST/$FILE ]]; then
		qsub -v R1=$SOURCE/$FILE $SCRIPT
	fi
done
