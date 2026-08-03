"""
Read all the files, and remove any taxa provided by the user.

This is relatively straightforward since the taxa files have the 
complete lineage for each taxon, so we can just check if any of the taxa
in the lineage match the taxa to be removed. If they do, then we skip
that taxon.

Our data has two columns - raw and normalised. We are not going to adjust
the normalised values, since we are only removing a small number of taxa, 
and we would need to know a priori which taxa to remove in order to adjust the normalised values.

Some of the (especially higher order) files will be the same!

Once this is run, we need to run

python ~/GitHubs/atavide_lite/summarise_taxonomy/scripts/join_taxonomies.py -t taxonomy -o taxonomy_summary_selected/ -n $SAMPLENAME

"""

import os
import sys
import argparse
import gzip
from atavide_lib import colours
__author__ = 'Rob Edwards'





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-d', '--directory', help='taxonomy directory with each sample', required=True)
    parser.add_argument('-o', '--output', help='output directory', required=True)
    parser.add_argument('-t', '--taxa', help='taxonomy to remove, e.g. f__Pseudomonadaceae', action='append', required=True)
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    for sample in os.listdir(args.directory):
        os.makedirs(os.path.join(args.output, sample), exist_ok=True)
        if args.verbose:
            print(f"{colours.OKBLUE}Processing sample {sample}{colours.ENDC}")
        for filename in os.listdir(os.path.join(args.directory, sample)):
            opener = open
            if filename.endswith('.gz'):
                opener = gzip.open
            with opener(os.path.join(args.directory, sample, filename), 'rt') as infile, opener(os.path.join(args.output, sample, filename), 'wt') as outfile:
                for line in infile:
                    if any(taxa in line for taxa in args.taxa):
                        continue
                    outfile.write(line)

print f"""
{colours.OKGREEN}Done!{colours.ENDC}
You can now run:
    python ~/GitHubs/atavide_lite/summarise_taxonomy/scripts/join_taxonomies.py -t {args.output} -o {args.output}_summary -n $SAMPLENAME
to create the summary files.
"""

