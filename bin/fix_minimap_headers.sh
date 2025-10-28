#!/bin/bash

# Fix the minimap2 headers issue
#
# If you have a large fasta file, minimap2 uses a multi-part index and does not (always) create the correct headers.
# The fix is a bit crap, you need to conver the bam to sam, and uncompress the fasta file, and then add the headers. 
# We do that using two temp files.
#
# Usage: fix_minimap_headers.sh <input.bam> <reference.fasta.gz> 

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.bam> <reference.fasta.gz>"
    exit 1
fi

input_bam="$1"
reference_fasta_gz="$2"

# uncompress the fasta file to a temp file
temp_fasta=$(mktemp --suffix=.fasta)
gunzip -c "$reference_fasta_gz" > "$temp_fasta"

# convert the bam to sam, and add the headers from the fasta file
temp_bam=$(mktemp --suffix=.bam)
samtools view -@16 $input_bam | samtools view -@16 -b -T "$temp_fasta" - > $temp_bam

# backup the original bam file
mv "$input_bam" "${input_bam}.bak"

# move the temp bam to the original bam file
mv "$temp_bam" "$input_bam"

# clean up the temp fasta file
rm "$temp_fasta"

echo "Fixed headers in $input_bam, original file backed up as ${input_bam}.bak"
