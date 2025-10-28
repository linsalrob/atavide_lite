#!/bin/bash
#SBATCH --job-name=vamb_mags
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -o slurm_output/vamb_slurm/%x-%j.out
#SBATCH -e slurm_output/vamb_slurm/%x-%j.err

set -euo pipefail

SRC="$HOME/GitHubs/atavide_lite/pawsey_shortread/"

VAMBUNSPLIT=$(sbatch --parsable \
	"$SRC/vamb_split.slurm" \
	vamb/contigs.fna.gz \
	vamb/vamb/vae_clusters_unsplit.tsv \
	20000 \
	vamb/vamb/clusters_unsplit)

echo "  -> UNSPLIT job ID: $VAMBUNSPLIT"

VAMBSPLIT=$(sbatch --parsable \
	"$SRC/vamb_split.slurm" \
	vamb/contigs.fna.gz \
	vamb/vamb/vae_clusters_split.tsv \
	20000 \
	vamb/vamb/clusters_split)

echo "  -> SPLIT job ID: $VAMBSPLIT"
