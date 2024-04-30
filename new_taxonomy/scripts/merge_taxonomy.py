"""
Merge the LCA information with the taxonkit taxonomy.
"""

import os
import sys
import gzip
import argparse
__author__ = 'Rob Edwards'


def parse_taxonomy(taxonfile, verbose=False):
    if verbose:
        print(f"Parsing taxonomy file {taxonfile}", file=sys.stderr)

    taxon = {}
    if taxonfile.endswith('.gz'):
        fin = gzip.open(taxonfile, 'rt')
    else:
        fin = open(taxonfile, 'r')

    for l in fin:
        p = l.rstrip().split("\t")
        if len(p) < 2:
            continue
        taxon[p[0]] = p[1]

    fin.close()

    return taxon


def write_lca(lcafile, taxa, outputfile, verbose=False):
    """
    Write the combined taxonomy out
    """

    if verbose:
        print(f"Writing {outputfile}", file=sys.stderr)

    if lcafile.endswith('.gz'):
        fin = gzip.open(lcafile, 'rt')
    else:
        fin = open(lcafile, 'r')

    if outputfile.endswith('.gz'):
        out = gzip.open(outputfile, 'wt')
    else:
        out = open(outputfile, 'w')
    
    for l in fin:
        p = l.strip().split("\t")
        if p[1] not in taxa:
            if verbose:
                print(f"No taxonomy for {p[1]}", file=sys.stderr)
            taxa[p[1]] = "k__;p__;c__;o__;f__;g__;s__";
        print("\t".join([l.rstrip(), taxa[p[1]]]), file=out)
    
    out.close()
    fin.close()

if __name__ == "__main2__":
    parser = argparse.ArgumentParser(description='Merge the taxonomy and the LCA data')
    parser.add_argument('-l', '--lca', help='LCA file from mmseqs easy-taxonomy', required=True)
    parser.add_argument('-t', '--taxonomy', help='taxonomy tuple file that has [taxid, taxonomy string] from taxonkit', required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()


    #taxa = parse_taxonomy(args.taxonomy, args.verbose)
    #write_lca(args.lca, taxa, args.output, args.verbose)


if __name__ == "__main__":
    taxa = parse_taxonomy(snakemake.input.tk)
    write_lca(snakemake.input.lca, taxa, snakemake.output.lcatax)


