#!/bin/bash
#SBATCH --job-name=vamb_split
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH -o vamb_split-%j.out
#SBATCH -e vamb_split-%j.err

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate bioinformatics

# use our local BGFS file system

if [[ $# -eq 0 ]]; then 
	echo "$0 <vamb contigs> <clusters file> <output directory>"
	exit;
fi

cp $1 $BGFS
CONTIGS=$(basename $1)
cp $2 $BGFS
CLUSTERS=$(basename $2)
ODIR=$(basename $3)
mkdir -p $BGFS/$ODIR

# This is the base command:
# python ~/GitHubs/atavide_lite/bin/vamb_create_fasta_clusters.py -f ../vamb/contigs.fna.gz -c vae_clusters_unsplit.tsv -o clusters_unsplit -v

echo "Running: $HOME/GitHubs/atavide_lite/bin/vamb_create_fasta_clusters.py -f $BGFS/$CONTIGS -c $BGFS/$CLUSTERS -o $BGFS/$ODIR" >&2

python $HOME/GitHubs/atavide_lite/bin/vamb_create_fasta_clusters.py -f $BGFS/$CONTIGS -c $BGFS/$CLUSTERS -o $BGFS/$ODIR

rsync -a $BGFS/$ODIR/ $3/
