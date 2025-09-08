#!/bin/bash
#SBATCH --job-name=vamb_mags
#SBATCH --time=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -o slurm_output/vamb_slurm/%x-%j.out
#SBATCH -e slurm_output/vamb_slurm/%x-%j.err

set -euo pipefail


E_BADARGS=65

if [[ "$#" -ne 1 ]]; then
	echo "$0 <file with [group, sample]>";
	exit $E_BADARGS;
fi

IFS=$'\n'
SRC="$HOME/GitHubs/atavide_lite/pawsey_shortread/"

for LN in $(grep -v Sample $1); do
	GRP=$(echo $LN | cut -d$'\t' -f 1)
	SMP=$(echo $LN | cut -d$'\t' -f 2)
	# check how many reads we have for $SMP. Should be only 1

	echo "Submitting UNSPLIT and SPLIT jobs for group: $GRP using $SRC/vamb_split.slurm"

	VAMBUNSPLIT=$(sbatch --parsable \
		"$SRC/vamb_split.slurm" \
		vamb_groups/$GRP/contigs.fna.gz \
		vamb_groups/$GRP/vamb/vae_clusters_unsplit.tsv \
		20000 \
		vamb_groups/$GRP/vamb/clusters_unsplit)

	echo "  -> UNSPLIT job ID: $VAMBUNSPLIT"

	VAMBSPLIT=$(sbatch --parsable \
		"$SRC/vamb_split.slurm" \
		vamb_groups/$GRP/contigs.fna.gz \
		vamb_groups/$GRP/vamb/vae_clusters_split.tsv \
		20000 \
		vamb_groups/$GRP/vamb/clusters_split)

	echo "  -> SPLIT job ID: $VAMBSPLIT"
done
