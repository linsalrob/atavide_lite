"""
Consume the multi-line json output and then summarise the results
"""

import os
import sys
import argparse
import json
__author__ = 'Rob Edwards'


def read_json(filename, verbose=False):
    if verbose:
        print(f"Reading {filename}", file=sys.stderr)
    with open(filename, 'r') as f:
        for l in f:
            try:
                js = json.loads(l)
            except json.JSONDecodeError as e:
                print(f"Error decoding JSON from line: {l.strip()}", file=sys.stderr)
                continue
            yield js

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-j', '--json', help='input json file. One record per line', required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    if not os.path.exists(args.json):
        print(f"Input file {args.json} does not exist", file=sys.stderr)
        sys.exit(1)

    with open(args.output, 'w') as out:
        print("Sample\t16S hits\tTotal sequences\t% 16S", file=out)
        for j in read_json(args.json, verbose=args.verbose):
            for sample in j:
                print(f"{sample}\t{j[sample]['QC-passed reads']['mapped']}\t{j[sample]['QC-passed reads']['total']}\t{j[sample]['QC-passed reads']['mapped %']}", file=out)

    


