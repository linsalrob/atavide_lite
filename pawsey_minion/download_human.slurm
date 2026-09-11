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

DEST=/scratch/$PAWSEY_PROJECT/$USER/Databases/human

mkdir -p $DEST
cd $DEST

if [[ -e GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz ]]; then
	echo "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz already exists. Nothing to do";
	touch GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz;
	exit 0
fi

curl -LO https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
echo "Downloaded human db to $DEST"
