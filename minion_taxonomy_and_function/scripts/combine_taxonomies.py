"""
concatentate all the ORF taxonomy files into 
a set of kingdom, phylum, etc files
"""

import os
import sys
import gzip
__author__ = 'Rob Edwards'



def concatenate_taxonomies(snakes):
    """
    Read each file and count the number of reads and the fraction for each type
    """

    # k__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__;s__
    kdm = {}
    phy = {}
    cls = {}
    odr = {}
    fam = {}
    gen = {} 
    spc = {}

    akdm = set()
    aphy = set()
    acls = set()
    aodr = set()
    afam = set()
    agen = set()
    aspc = set()

    counts = {}

    for t in snakes.input.taxfiles:
        sample = os.path.basename(t).replace(".contigs.tsv.gz", "")
        opener = open
        if t.endswith('.gz'):
            opener = gzip.open
        kdm[sample] = {}
        phy[sample] = {}
        cls[sample] = {}
        odr[sample] = {}
        fam[sample] = {}
        gen[sample] = {}
        spc[sample] = {}

        counts[sample] = 0
        with opener(t, 'rt') as f:
            for l in f:
                counts[sample] += 1
                p = l.strip().split("\t")
                if len(p) != 8:
                    print(f"ERROR: if file {t} this line doesn't have the right columns\n{l}", file=sys.stderr)
                    exit(1)
                kdm[sample][p[1]] =  kdm[sample].get(p[1], 0)+1; akdm.add(p[1])
                phy[sample][p[2]] =  phy[sample].get(p[2], 0)+1; aphy.add(p[2])
                cls[sample][p[3]] =  cls[sample].get(p[3], 0)+1; acls.add(p[3])
                odr[sample][p[4]] =  odr[sample].get(p[4], 0)+1; aodr.add(p[4])
                fam[sample][p[5]] =  fam[sample].get(p[5], 0)+1; afam.add(p[5])
                gen[sample][p[6]] =  gen[sample].get(p[6], 0)+1; agen.add(p[6])
                spc[sample][p[7]] =  spc[sample].get(p[7], 0)+1; aspc.add(p[7])

    samples = sorted(kdm.keys())

    datasets = (
        ("kingdom", kdm, akdm),
        ("phylum", phy, aphy),
        ("class", cls, acls),
        ("order", odr, aodr),
        ("family", fam, afam),
        ("genus", gen, agen),
        ("species", spc, aspc)
    )

    os.makedirs(snakes.output.outdir, exist_ok=True)
    for tple in datasets:
        with gzip.open(os.path.join(snakes.output.outdir, f"{tple[0]}.tsv.gz"), 'wt') as out:
            headers = list(sorted(tple[2]))
            c = tple[1]
            print("\t".join(map(str, ["Sample"] + headers)), file=out)
            for s in samples:
                print(s, file=out, end="")
                for h in headers:
                    if h in c[s]:
                        norm = (c[s][h] / counts[s]) * 1e6
                        print(f"\t{norm}", end="", file=out)
                    else:
                        print("\t0", end="", file=out)
                print(file=out)




concatenate_taxonomies(snakemake)
