#!/bin/bash
#SBATCH --job-name=human_download
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1800M
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err

# NOTE: Setonix enforces MaxMemPerCPU (work=1840M, highmem=7900M). A bare --mem
# larger than cpus-per-task * MaxMemPerCPU does not fail - Slurm silently raises the
# CPU count to cover it, so the job queues for far more CPUs than it uses. Request
# memory with --mem-per-cpu so the two can never disagree.

set -euo pipefail

DEST=/scratch/$PAWSEY_PROJECT/$USER/Databases/human

mkdir -p $DEST
cd $DEST

if gzip -t GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz 2>/dev/null; then
	echo "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz already exists. Nothing to do";
	touch GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz;
	exit 0
fi

# NCBI moved this file; the old GCF_000001405.40 .../GRCh38_major_release_seqs_for_alignment_pipelines/
# path now returns an HTML page with HTTP 404, which a plain "curl -LO" saved as the .fna.gz.
# -f makes curl fail on an HTTP error, and gzip -t rejects anything that is not a whole gzip file.
F=GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
curl -fL -o "$F.part" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/$F"
gzip -t "$F.part"
mv "$F.part" "$F"
echo "Downloaded human db to $DEST"
