"""
Summarise ORF-based taxonomy, including based on the location in the
genome
"""

import os
import sys
import argparse
import gzip
import re


__author__ = 'Rob Edwards'


def parse_taxonomy(taxfile, outputfile):
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
            kng, phy, cls, orr, fam, gen, sp = p[4].split(";")

            ctgre = re.compile(r'^(.*)-(orf\d+)')
            m = ctgre.match(p[0])
            if not m:
                print(f"Error: no group from {p[0]} in {taxfile}", file=sys.stderr)
                continue
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


    
    # for each contig we have the best hit, so we can tax the max for tax!
    print("\t".join(['contig', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']))



    openwriter = open
    if outputfile.endswith('.gz'):
        openwriter = gzip.open

    with openwriter(outputfile, 'wt') as out:
        for contig in contig_tax.keys():
            n+=1
            print(contig, end="", file=out)
            for tx in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
                contigmx = max(contig_tax[contig][tx], key=contig_tax[contig][tx].get)
                print(f"\t{contigmx}", end="", file=out)
                tax[tx][contigmx] = tax[tx].get(contigmx, 0)+1
            print(file=out)

    # now we have the counts.
    print("Most likely taxonomy overall", file=sys.stderr)
    for tx in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
        mx = max(tax[tx], key=tax[tx].get)
        perc = tax[tx][mx]/n
        print(f"{tx}\t{mx}\t{perc}", file=sys.stderr)



    

parse_taxonomy(snakemake.input.lcatax, snakemake.output.consum)



