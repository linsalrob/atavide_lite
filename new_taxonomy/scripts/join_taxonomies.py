"""
Merge all the taxonomies into a set of files. I can not figure out how to
do this with snakemake checkpoints so I gave up!
"""
import argparse
import os
import sys
import gzip
import re

__author__ = 'Rob Edwards'

def join(files, raw_output, norm_output):
    raw = {}
    norm = {}
    sample_names = {}

    alltax = set()
    for f in files:
        raw[f] = {}
        norm[f] = {}

        opener = gzip.open if f.endswith('.gz') else open
        with opener(f, 'rt') as fn:
            for l in fn:
                p = l.strip().split("\t")
                raw[f][p[2]] = p[3]
                norm[f][p[2]] = p[4]
                sample_names[f] = p[1]
                alltax.add(p[2])

    with gzip.open(raw_output, 'wt') as out:
        print("\t".join(["taxonomy"] + files), file=out)
        for t in sorted(alltax):
            print(t, file=out, end="")
            for f in files:
                if t not in raw[f]:
                    print("\t0", file=out, end="")
                else:
                    print(f"\t{raw[f][t]}", file=out, end="")
            print(file=out)

    with gzip.open(norm_output, 'wt') as out:

        print("\t".join(["taxonomy"] + [sample_names[f] for f in files]), file=out)
        for t in sorted(alltax):
            print(t, file=out, end="")
            for f in files:
                if t not in norm[f]:
                    print("\t0", file=out, end="")
                else:
                    print(f"\t{norm[f][t]}", file=out, end="")
            print(file=out)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Join the taxonomies into separate files')
    parser.add_argument('-t', '--taxonomy', help='taxonomy directory with all the samples', required=True)
    parser.add_argument('-n', '--name', help='sample name to add to the output files', default='atavide')
    parser.add_argument('-o', '--output', help='output directory', required=True)
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()
    taxa = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    os.makedirs(args.output, exist_ok=True)
    samplename = re.sub(r'\W', '', args.name)

    for t in taxa:
        inputs = []
        for sd in os.listdir(args.taxonomy):
            if os.path.exists(os.path.join(args.taxonomy, sd, f"{t}.tsv.gz")):
                inputs.append(os.path.join(args.taxonomy, sd, f"{t}.tsv.gz"))
        if args.verbose:
            print(f"Writing {t} from {inputs}", file=sys.stderr)
        join(inputs, os.path.join(args.output, f"{samplename}_{t}.raw.tsv.gz"),  os.path.join(args.output, f"{samplename}_{t}.norm.tsv.gz"))
