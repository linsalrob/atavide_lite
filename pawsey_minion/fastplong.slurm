#!/bin/bash
# Long-read QC with fastplong. Submit as a Slurm array over reads.txt.
# Create slurm_output/fastplong_slurm before submitting, since Slurm opens
# stdout/stderr before this script begins.
#
# Example:
#   NUM_READS=$(wc -l < reads.txt)
#   mkdir -p slurm_output/fastplong_slurm
#   sbatch --array=1-${NUM_READS}%4 pawsey_minion/fastplong.slurm
#
# Output FASTQs are intentionally written to fastq_fastplong/ rather than
# fastq_fastp/ so this revised QC run cannot overwrite the original results.

#SBATCH --job-name=fastplong
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm_output/fastplong_slurm/fastplong-%A_%a.out
#SBATCH --error=slurm_output/fastplong_slurm/fastplong-%A_%a.err

set -euo pipefail

if [[ ! -f DEFINITIONS.sh ]]; then
    echo "ERROR: run from the analysis directory containing DEFINITIONS.sh." >&2
    exit 2
fi
if [[ ! -f reads.txt ]]; then
    echo "ERROR: reads.txt is required for the Slurm array lookup." >&2
    exit 2
fi
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: submit this script as a Slurm array." >&2
    exit 2
fi

source DEFINITIONS.sh

READS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" reads.txt)
if [[ -z "$READS" ]]; then
    echo "ERROR: no reads.txt entry for array task ${SLURM_ARRAY_TASK_ID}." >&2
    exit 2
fi
if [[ ! -f "fastq/$READS" ]]; then
    echo "ERROR: input FASTQ not found: fastq/$READS" >&2
    exit 2
fi

BASE=${READS%"$FILEEND"}
if [[ "$BASE" == "$READS" ]]; then
    echo "ERROR: $READS does not end with FILEEND=$FILEEND." >&2
    exit 2
fi

CONDA_PREFIX=/scratch/${PAWSEY_PROJECT}/${USER}/software/miniconda3/atavide_lite
FASTPLONG="$CONDA_PREFIX/bin/fastplong"
ADAPTERS="$HOME/GitHubs/atavide_lite/adapters/ONT_RPB114.24.fna"

if [[ ! -x "$FASTPLONG" ]]; then
    echo "ERROR: fastplong executable not found: $FASTPLONG" >&2
    exit 2
fi
if [[ ! -f "$ADAPTERS" ]]; then
    echo "ERROR: ONT adapter FASTA not found: $ADAPTERS" >&2
    exit 2
fi

mkdir -p fastq_fastplong fastplong_output slurm_output/fastplong_slurm

OUTPUT="fastq_fastplong/$READS"
if [[ -e "$OUTPUT" ]]; then
    echo "ERROR: refusing to overwrite existing output: $OUTPUT" >&2
    exit 2
fi

# Retain usable long-read sequence while avoiding aggressive global read removal:
# -Q disables the global quality filter; --break removes low-quality sections
# within reads. Strict adapter matching avoids spurious adapter trimming.
"$FASTPLONG" \
    -i "fastq/$READS" \
    -o "$OUTPUT" \
    --adapter_fasta "$ADAPTERS" \
    --distance_threshold 0.1 \
    --trimming_extension 0 \
    -Q \
    --break \
    --break_window_size 100 \
    --break_mean_quality 10 \
    --n_base_limit 2 \
    -l 20 \
    -w "$SLURM_CPUS_PER_TASK" \
    -j "fastplong_output/${BASE}.json" \
    -h "fastplong_output/${BASE}.html"
