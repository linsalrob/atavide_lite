"""
DEPRECATED: This script depends on the VAMB Python module and is subject to
breaking changes between VAMB versions.

Please use vamb_create_fasta_clusters.py instead, which is VAMB-version
independent and does not require importing the vamb module.

Usage:
    python bin/vamb_create_fasta_clusters.py -f FASTA -c CLUSTERS -o OUTDIR [-m MINSIZE]

This script is kept for backwards compatibility but may not work with all
VAMB versions. See docs/compat.md for more information about VAMB compatibility.
"""

import sys
import argparse
import pathlib

# Check VAMB version before importing
try:
    import vamb
    try:
        VAMB_VERSION = vamb.__version__
    except AttributeError:
        VAMB_VERSION = "unknown"
        print("Warning: Could not determine VAMB version", file=sys.stderr)
        print("This script is tested with VAMB 4.0.0+", file=sys.stderr)
except ImportError:
    print("Error: VAMB module not found", file=sys.stderr)
    print("", file=sys.stderr)
    print("RECOMMENDED: Use vamb_create_fasta_clusters.py instead,", file=sys.stderr)
    print("which does not require the VAMB module.", file=sys.stderr)
    print("", file=sys.stderr)
    print("If you want to use this script, install VAMB:", file=sys.stderr)
    print("  conda install vamb=4.1.0", file=sys.stderr)
    sys.exit(1)

# Warn if not tested version
if VAMB_VERSION != "unknown":
    major_version = int(VAMB_VERSION.split('.')[0]) if '.' in VAMB_VERSION else 0
    if major_version < 4:
        print(f"Warning: VAMB {VAMB_VERSION} detected. This script is tested with VAMB 4.x", file=sys.stderr)
        print("Consider using vamb_create_fasta_clusters.py for better compatibility", file=sys.stderr)

parser = argparse.ArgumentParser(
    description="""Command-line bin creator.
Will read the entire content of the FASTA file into memory - beware.""",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    add_help=False,
)

parser.add_argument("fastapath", help="Path to FASTA file")
parser.add_argument("clusterspath", help="Path to clusters.tsv")
parser.add_argument("minsize", help="Minimum size of bin in bp", type=int, default=0)
parser.add_argument("outdir", help="Directory to create")

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

clusters = {
    cluster: contigs
    for (cluster, contigs) in clusters.items()
    if sum(lens[c] for c in contigs) >= args.minsize
}

if not clusters:
    print(f"No bins passed the minimum size threshold of {minsize}", file=sys.stderr)
    sys.exit(0)

print(f"Writing {len(clusters)} bins to {args.outdir}", file=sys.stderr)
with vamb.vambtools.Reader(args.fastapath) as file:
    vamb.vambtools.write_bins(pathlib.Path(args.outdir), clusters, file, maxbins=None, compress_output=True)
