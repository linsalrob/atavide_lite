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


join(snakemake.input.files, snakemake.output.raw_output, snakemake.output.norm_output)
