#!/usr/bin/env python3
"""Remove FASTQ records named in a spike-in alignment list."""

import argparse
import gzip


def opener(path, mode):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


parser = argparse.ArgumentParser()
parser.add_argument("--ids", required=True)
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--stats", required=True)
args = parser.parse_args()

with open(args.ids, encoding="utf-8") as handle:
    spike_ids = {line.rstrip("\n") for line in handle if line.strip()}

total = removed = 0
with opener(args.input, "rt") as source, opener(args.output, "wt") as dest:
    while True:
        record = [source.readline() for _ in range(4)]
        if not record[0]:
            if any(record[1:]):
                raise ValueError("truncated FASTQ record at end of input")
            break
        if not all(record) or not record[0].startswith("@"):
            raise ValueError("invalid or truncated FASTQ record")
        total += 1
        read_id = record[0][1:].split(maxsplit=1)[0]
        if read_id in spike_ids:
            removed += 1
        else:
            dest.writelines(record)

with open(args.stats, "w", encoding="utf-8") as output:
    output.write("input_reads\tremoved_spike_in_reads\tretained_reads\tunique_spike_in_read_ids\n")
    output.write(f"{total}\t{removed}\t{total - removed}\t{len(spike_ids)}\n")
