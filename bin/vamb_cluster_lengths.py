import sys
import argparse
import vamb
import pathlib

parser = argparse.ArgumentParser(
    description="""Command-line bin creator.
Will read the entire content of the FASTA file into memory - beware.""",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    add_help=False,
)

parser.add_argument("fastapath", help="Path to FASTA file")
parser.add_argument("clusterspath", help="Path to clusters.tsv")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()

args = parser.parse_args()

# Read in FASTA files only to get its length. This way, we can avoid storing
# in memory contigs for sequences that will never get output anyway
lens: dict[str, int] = dict()
with vamb.vambtools.Reader(args.fastapath) as file:
    for record in vamb.vambtools.byte_iterfasta(file, args.fastapath):
        lens[record.identifier] = len(record)

with open(args.clusterspath) as file:
    clusters = vamb.vambtools.read_clusters(file)

cluster_lengths = {}
for (cluster, contigs) in clusters.items():
    cluster_lengths[cluster] = sum(lens[c] for c in contigs)

cluster_lengths_count = {}
for length in cluster_lengths.values():
    if length not in cluster_lengths_count:
        cluster_lengths_count[length] = 0
    cluster_lengths_count[length] += 1

for c in sorted(cluster_lengths_count.keys()):
    print(f"{c}\t{cluster_lengths_count[c]}")



