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


def parse_taxonomy(taxfile, output_directory):
    """
    Parse the taxonomy file
    """

    opener = open

    if ".gz" in taxfile:  # gzipped version of the file, we need gzip.open
        opener = gzip.open

    tax = {
        'kingdom'  : {},
        'phylum'  : {},
        'class'  : {},
        'order'  : {},
        'family'  : {},
        'genus'  : {},
        'species'  : {}
    }

    readn = 0
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
            kng, phy, cls, orr, fam, gen, sp = p[8].split(";")


            tax['kingdom'][kng] = tax['kingdom'].get(kng, 0) + 1
            tax['phylum'][f'{kng};{phy}'] = tax['phylum'].get(f'{kng};{phy}', 0) + 1
            tax['class'][f'{kng};{phy};{cls}'] = tax['class'].get(f'{kng};{phy};{cls}', 0) + 1
            tax['order'][f'{kng};{phy};{cls};{orr}'] = tax['order'].get(f'{kng};{phy};{cls};{orr}', 0) + 1
            tax['family'][f'{kng};{phy};{cls};{orr};{fam}'] = tax['family'].get(f'{kng};{phy};{cls};{orr};{fam}', 0) + 1
            tax['genus'][f'{kng};{phy};{cls};{orr};{fam};{gen}'] = tax['genus'].get(f'{kng};{phy};{cls};{orr};{fam};{gen}', 0) + 1
            tax['species'][f'{kng};{phy};{cls};{orr};{fam};{gen};{sp}'] = tax['species'].get(f'{kng};{phy};{cls};{orr};{fam};{gen};{sp}', 0) + 1
            readn += 1



    os.makedirs(output_directory, exist_ok=True)
    sampleid = os.path.basename(output_directory)
    for tx in tax.keys():
        with gzip.open(os.path.join(output_directory, f"{tx}.tsv.gz"), 'wt') as out:
            for nm in tax[tx]:
                normed = tax[tx][nm]/readn * 1000000
                print(f"{tx}\t{sampleid}\t{nm}\t{tax[tx][nm]}\t{normed}", file=out)
    


    

parse_taxonomy(snakemake.input.lcatax, snakemake.output.directory)



