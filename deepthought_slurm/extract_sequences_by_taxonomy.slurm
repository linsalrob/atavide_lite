#!/bin/bash
#SBATCH --time=10-0
#SBATCH --job-name=extract_seqs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -o slurm_output/extract_sequences-%j.out
#SBATCH -e slurm_output/extract_sequences-%j.err

if [[ $# -ne 2 ]]; then 
	echo "$0 <taxonomy> <outputdir>"
	exit; 
fi

echo "Writing outputs to $BGFS/$2"

python /home/edwa0468/GitHubs/atavide_lite/bin/extract_reads_by_taxonomy.py -t $1 -o $BGFS/$2 -v
cp -r $BGFS/$2 $2

