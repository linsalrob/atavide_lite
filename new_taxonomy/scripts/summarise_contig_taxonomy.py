"""
Summarise ORF-based taxonomy, including based on the location in the
genome
"""

import os
import sys
import argparse
import gzip
import re

from roblib import colours, stream_fasta

__author__ = 'Rob Edwards'


def orf_locations(orf_file, verbose=False):
    """
    Read the ORF locations and store the midpoint
    """

    if verbose:
        print(f"{colors.BLUE}Reading {orf_file}{colours.ENDC}", 
              file=sys.stderr)


    opener = open

    if ".gz" in orf_file:  # gzipped version of the file, we need gzip.open
        opener = gzip.open

    locs = re.compile(r'>(\w+)\s+\[[\w\-]+\s+frame\s+[\+\-\d+]+\s(\d+)\s(\d+)\]')

    position = {}

    with opener(orf_file, "rt") as f:
        # >orf1 [4a901472-3d00-4da0-b0e7-9169d21efcfc frame +1 58 225]
        for l in f:
            m = locs.match(l)
            (orf, beg, end) = m.groups()
            position[orf] = (int(beg)+int(end))/2
    return position


def parse_taxonomy(taxfile, verbose=False):
    """
    Parse the taxonomy file
    """

    opener = open

    if ".gz" in taxfile:  # gzipped version of the file, we need gzip.open
        opener = gzip.open

    tax = {
        'kingdom' : {},
        'phylum' : {},
        'class' : {},
        'order' : {},
        'family' : {},
        'genus' : {},
        'species' : {}
    }
    contig_tax = {}
    contign = {}
    n = 0
    with opener(taxfile, "rt") as f:
        for l in f:
            p = l.strip().split("\t")
            if p[1] == '0' or p[1] == '131567' or p[3] == "unclassified":
                continue
            if int(p[1]) < 3:
                continue
            # We have this data
            # k__;p__;c__;o__;f__;g__;s__
            #
            # but we need to keep the components together for the counts
            # because it won't make any sense if we merge 
            # just based on the most abundant
            #
            #
            if verbose:
                print(l, file=sys.stderr, end="")
            kng, phy, cls, orr, fam, gen, sp = p[4].split(";")

            ctgre = re.compile(r'^(.*)-(orf\d+)')
            m = ctgre.match(p[0])
            contig = m.groups()[0]

            if contig not in contig_tax:
                contig_tax[contig] = {
                    'kingdom' : {},
                    'phylum' : {},
                    'class' : {},
                    'order' : {},
                    'family' : {},
                    'genus' : {},
                    'species' : {}
                }

            contig_tax[contig]['kingdom'][kng] = contig_tax[contig]['kingdom'].get(kng, 0) + 1
            contig_tax[contig]['phylum'][f'{kng};{phy}'] = contig_tax[contig]['phylum'].get(f'{kng};{phy}', 0) + 1
            contig_tax[contig]['class'][f'{kng};{phy};{cls}'] = contig_tax[contig]['class'].get(f'{kng};{phy};{cls}', 0) + 1
            contig_tax[contig]['order'][f'{kng};{phy};{cls};{orr}'] = contig_tax[contig]['order'].get(f'{kng};{phy};{cls};{orr}', 0) + 1
            contig_tax[contig]['family'][f'{kng};{phy};{cls};{orr};{fam}'] = contig_tax[contig]['family'].get(f'{kng};{phy};{cls};{orr};{fam}', 0) + 1
            contig_tax[contig]['genus'][f'{kng};{phy};{cls};{orr};{fam};{gen}'] = contig_tax[contig]['genus'].get(f'{kng};{phy};{cls};{orr};{fam};{gen}', 0) + 1
            contig_tax[contig]['species'][f'{kng};{phy};{cls};{orr};{fam};{gen};{sp}'] = contig_tax[contig]['species'].get(f'{kng};{phy};{cls};{orr};{fam};{gen};{sp}', 0) + 1
            contign[contig] = contign.get(contig, 0) + 1

            """
            
            
            # now we reconstruct these at all levels and count them
            tax['kingdom'][kng] = tax['kingdom'].get(kng, 0) + 1
            tax['phylum'][f'{kng};{phy}'] = tax['phylum'].get(f'{kng};{phy}', 0) + 1
            tax['class'][f'{kng};{phy};{cls}'] = tax['class'].get(f'{kng};{phy};{cls}', 0) + 1
            tax['order'][f'{kng};{phy};{cls};{orr}'] = tax['order'].get(f'{kng};{phy};{cls};{orr}', 0) + 1
            tax['family'][f'{kng};{phy};{cls};{orr};{fam}'] = tax['family'].get(f'{kng};{phy};{cls};{orr};{fam}', 0) + 1
            tax['genus'][f'{kng};{phy};{cls};{orr};{fam};{gen}'] = tax['genus'].get(f'{kng};{phy};{cls};{orr};{fam};{gen}', 0) + 1
            tax['species'][f'{kng};{phy};{cls};{orr};{fam};{gen};{sp}'] = tax['species'].get(f'{kng};{phy};{cls};{orr};{fam};{gen};{sp}', 0) + 1
            
            orf_tax[p[0]] = [kng, phy, cls, orr, fam, gen, sp]
            """

    
    # for each contig we have the best hit, so we can tax the max for tax!
    print("\t".join(['contig', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']))

    for contig in contig_tax.keys():
        n+=1
        print(contig, end="")
        for tx in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
            contigmx = max(contig_tax[contig][tx], key=contig_tax[contig][tx].get)
            print(f"\t{contigmx}", end="")
            tax[tx][contigmx] = tax[tx].get(contigmx, 0)+1
        print()

    # now we have the counts.
    print("Most likely taxonomy overall", file=sys.stderr)
    for tx in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
        mx = max(tax[tx], key=tax[tx].get)
        perc = tax[tx][mx]/n
        print(f"{tx}\t{mx}\t{perc}", file=sys.stderr)



    





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-o', help='orf fasta file that includes the gene locations', required=True)
    parser.add_argument('-t', help='lca taxonomy file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    
    parse_taxonomy(args.t, args.v)


