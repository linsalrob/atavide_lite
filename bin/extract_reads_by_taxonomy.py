"""
We are going to read the mmseqs taxonomy files, look for a particular taxonomy, and then extract the reads that match that taxonomy.
"""
import gzip
import os
import sys
import argparse

__author__ = 'Rob Edwards'

from atavide_lib import colors, read_definitions

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-t', '--taxonomy', help='taxonomy to extract', required=True)
    parser.add_argument('-o', '--output', help='output directory', required=True)
    parser.add_argument('-m', '--mmseqs', help='mmseqs directory', default='mmseqs')
    parser.add_argument('-f', '--fasta', help='fasta directory', default='fasta')
    parser.add_argument('-d', '--definitions', help='definitions file', default='DEFINITIONS.sh')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # read the definitions file
    if not os.path.exists('DEFINITIONS.sh'):
        print(f"{colors.RED}Could not find DEFINITIONS.sh{colors.ENDC}")
        sys.exit(1)

    defs = read_definitions(args.definitions)

    # read each of the files in the mmseqs directory
    for sample in os.listdir(args.mmseqs):
        if not os.path.exists(os.path.join(args.mmseqs, sample, f"{sample}_lca_taxonomy.tsv.gz")):
            print(f"{colors.RED}Could not find {sample}_lca_taxonomy.tsv.gz{colors.ENDC}")
            continue
        if args.verbose:
            print(f"{colors.GREEN}Reading {sample}_lca_taxonomy.tsv.gz{colors.ENDC}", file=sys.stderr)
        wanted = set()
        with gzip.open(os.path.join(args.mmseqs, sample, f"{sample}_lca_taxonomy.tsv.gz"), 'rt') as f:
            for l in f:
                p = l.strip().split("\t")
                if args.taxonomy in p[8]:
                    wanted.add(p[0])
        # now read the fasta file from the fasta dir and write out the wanted reads
        if args.verbose:
            print(f"{colors.BLUE}\tWriting {len(wanted)} sequences to the fasta file", file=sys.stderr)
        with gzip.open(os.path.join(args.fasta, f"{sample}{defs['FAFILEEND']}"), 'rt') as f, open(os.path.join(args.output, f"{sample}_{args.taxonomy}.fasta"), 'w') as out:
            write = False
            written = 0
            for l in f:
                if l.startswith(">"):
                    write = False
                    p = l[1:].strip().split()
                    if p[0] in wanted:
                        written += 1
                        write = True
                if write:
                    out.write(l)
        if args.verbose:
            print(f"{colors.PINK}Wrote {written} sequences{colors.ENDC}", file=sys.stderr)
