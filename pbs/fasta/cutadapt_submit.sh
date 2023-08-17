#!/bin/bash


if [[ $# != 1 ]]; then
	echo "$0 <directory of (split) fasta files>" >&2;
	exit 1;
fi

for R1 in $(find $1 -name \*_R1*); do 
	R2=${R1/_R1/_R2}
	if [[ -e $R2 ]]; then
		if [[ ! -e cutadapt/$R1 ]] && [[ ! -e cutadapt/$R2 ]]; then
			qsub -v R1=$R1,R2=$R2 ~/GitHubs/atavide_lite/pbs/fasta/cutadapt.pbs
		elif [[ -e  cutadapt/$R1 ]]; then 
			qsub -v R2=$R2 ~/GitHubs/atavide_lite/pbs/fasta/cutadapt.pbs
		elif [[ -e cutadapt/$R2 ]]; then
			qsub -v R1=$R1 ~/GitHubs/atavide_lite/pbs/fasta/cutadapt.pbs
		fi
	else
		if [[ ! -e  cutadapt/$R1 ]]; then
			qsub -v R1=$R1 ~/GitHubs/atavide_lite/pbs/fasta/cutadapt.pbs
		fi
	fi
done
