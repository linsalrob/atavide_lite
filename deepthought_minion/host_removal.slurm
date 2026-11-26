#!/bin/bash
#SBATCH --job-name=HostRemoval
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH -o slurm_output/host_slurm/host-%A_%a.out
#SBATCH -e slurm_output/host_slurm/host-%A_%a.err

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate atavide_lite


if [[ ! -e reads.txt ]]; then
	echo "Please make a file with the reads using this command:" >&2
	echo "find fastq -type f -printf "%f\n" > reads.txt" >&2;
	exit 2;
fi

if [[ ! -e DEFINITIONS.sh ]]; then
	echo "Please create a DEFINITIONS.sh file with SOURCE, FILEEND, HOSTREMOVED, HOSTFILE, etc" >&2
	exit 2;
fi

source DEFINITIONS.sh


QC=${QC_DIR:-fastq_fastp}
THREADS=${SLURM_CPUS_PER_TASK:-16}
READS=$(head -n $SLURM_ARRAY_TASK_ID reads.txt | tail -n 1)

echo "Processing $READS" >&2;

if [[ ! -f "$QC/$READS" ]]; then
	echo "ERROR: QC FASTQ not found: $QC/$READS" >&2;
	exit 2;
fi

# When a spike-in reference is configured (SPIKE_IN_SEQUENCE, and optionally SPIKE_IN,
# in DEFINITIONS.sh), keep the ordinary host-removed reads as an intermediate and make
# HOSTREMOVED the spike-free downstream directory, as pawsey_minion/host_removal.slurm
# does. Without a spike-in, this is the former workflow.
PRE_SPIKE_HOSTREMOVED=$HOSTREMOVED
if [[ -n "${SPIKE_IN_SEQUENCE:-}" ]]; then
	PRE_SPIKE_HOSTREMOVED="${HOSTREMOVED}_before_${SPIKE_IN:-spike_in}"
fi
mkdir -p $BGFS/$HOST $BGFS/$PRE_SPIKE_HOSTREMOVED

O=${READS/$FILEEND/.bam}
if [[ $O == $READS ]]; then
	echo "ERROR: OUTPUT $O and INPUT $READS are the same. Are you sure the files end $FILEEND?" >&2;
	exit 1;
fi

TMPDIR=$(mktemp -d -p $BGFS)
echo -e "minimap2 -t $THREADS --secondary=no --split-prefix=$TMPDIR/tmp$$ -a -x map-ont $HOSTFILE $QC/$READS | samtools view -bh | samtools sort -o $BGFS/$HOST/$O -\n" >&2;
minimap2 -t $THREADS --secondary=no --split-prefix=$TMPDIR/tmp$$ -a -x map-ont $HOSTFILE "$QC/$READS" | samtools view -bh | samtools sort -o $BGFS/$HOST/$O -
samtools index $BGFS/$HOST/$O


# Reads that match the host
samtools fastq -F 3588 $BGFS/$HOST/$O | gzip -c > $BGFS/$HOST/$READS


# Reads that do not match the host

samtools fastq -F 3584 -f 4 $BGFS/$HOST/$O | gzip -c > $BGFS/$PRE_SPIKE_HOSTREMOVED/$READS

rsync -a $BGFS/$PRE_SPIKE_HOSTREMOVED/ $PRE_SPIKE_HOSTREMOVED/
rsync -a $BGFS/$HOST/ $HOST/

if [[ -n "${SPIKE_IN_SEQUENCE:-}" ]]; then
	if [[ ! -s "$SPIKE_IN_SEQUENCE" ]]; then
		echo "ERROR: SPIKE_IN_SEQUENCE is missing or empty: $SPIKE_IN_SEQUENCE" >&2
		exit 2
	fi
	mkdir -p "$HOSTREMOVED" spike_in/bam spike_in/removal_stats
	SPIKE_BAM="spike_in/bam/${READS/$FILEEND/}.${SPIKE_IN:-spike_in}.bam"
	SPIKE_STATS="spike_in/removal_stats/${READS/$FILEEND}.tsv"
	if [[ -e "$HOSTREMOVED/$READS" || -e "$SPIKE_BAM" || -e "$SPIKE_STATS" ]]; then
		echo "ERROR: refusing to overwrite existing spike-in outputs for $READS" >&2
		exit 2
	fi
	SPIKE_TMP=$(mktemp -d -p "$BGFS")
	# Keep primary, confidently mapped (MAPQ >=20) spike-in reads only. This
	# avoids subtracting incidental low-confidence ONT alignments.
	minimap2 -t "$THREADS" --secondary=no -a -x map-ont "$SPIKE_IN_SEQUENCE" "$PRE_SPIKE_HOSTREMOVED/$READS" \
		| samtools view -@ "$THREADS" -bh -F 2308 -q 20 - \
		| samtools sort -@ "$THREADS" -o "$SPIKE_BAM" -
	samtools index "$SPIKE_BAM"
	samtools view "$SPIKE_BAM" | cut -f1 | sort -u > "$SPIKE_TMP/spike_read_ids.txt"
	python "$HOME/GitHubs/atavide_lite/bin/remove_spikein_reads.py" \
		--ids "$SPIKE_TMP/spike_read_ids.txt" \
		--input "$PRE_SPIKE_HOSTREMOVED/$READS" \
		--output "$HOSTREMOVED/$READS" \
		--stats "$SPIKE_STATS"
fi
