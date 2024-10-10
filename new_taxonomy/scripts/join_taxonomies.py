"""
Merge all the taxonomies into a set of files. I can not figure out how to
do this with snakemake checkpoints so I gave up!
"""
import argparse
import os
import sys
import gzip

__author__ = 'Rob Edwards'

def join(files, raw_output, norm_output):
    raw = {}
    norm = {}
    alltax = set()
    for f in files:
        raw[f] = {}
        norm[f] = {}
        opener = gzip.open if f.endswith('.gz') else open
        with opener(f, 'r') as fn:
            for l in fn:
                p = l.strip().split("\t")
                raw[f][p[0]] = p[1]
                norm[f][p[0]] = p[2]
                alltax.add(p[0])

    with gzip.open(raw_output, 'wt') as out:
        print("\t".join(["taxonomy"] + files), file=out)
        for t in sorted(alltax):
            print(t, file=out, end="")
            for f in files:
                if t not in raw[f]:
                    print("\t0", file=out, end="")
                else:
                    print(f"\t{raw[f][t]}", file=out, end="")

    with gzip.open(norm_output, 'wt') as out:
        print("\t".join(["taxonomy"] + files), file=out)
        for t in sorted(alltax):
            print(t, file=out, end="")
            for f in files:
                if t not in norm[f]:
                    print("\t0", file=out, end="")
                else:
                    print(f"\t{norm[f][t]}", file=out, end="")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Join the taxonomies into separate files')
    parser.add_argument('-t', '--taxonomy', help='taxonomy directory with all the samples', required=True)
    parser.add_argument('-o', '--output', help='output directory', required=True)
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()
    taxa = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    os.makedirs(args.output, exist_ok=True)

    for t in taxa:
        inputs = []
        for sd in os.listdir(args.taxonomy):
            if os.path.exists(os.path.join(args.taxonomy, sd, f"{t}.tsv.gz")):
                inputs.append(os.path.join(args.taxonomy, sd, f"{t}.tsv.gz"))
        if args.verbose:
            print(f"Writing {t} from {inputs}", file=sys.stderr)
        join(inputs, os.path.join(args.output, f"{t}.raw.tsv.gz"),  os.path.join(args.output, f"{t}.norm.tsv.gz"))