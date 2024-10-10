"""

"""

import os
import sys
import argparse
import gzip

__author__ = 'Rob Edwards'



def count_ss(ssfiles, ssdir):
    """
    Count the subsystems in ssfiles
    """

    read_count = {} # total of all reads
    ss_read_count = {} # total of only those reads in ss
    cls = {}
    lvl1 = {}
    lvl2 = {}
    ss = {}
    ssall = {}
    # note we hash "backwards" with the primary key the ss name and then the filename
    # so we don't need a separate list of ss names

    for ssf in ssfiles:
        opener = open
        if ssf.endswith('.gz'):
            opener = gzip.open
        sid = os.path.basename(ssf).replace("_tophit_report_subsystems.gz", "")

        read_count[sid] = 0
        ss_read_count[sid] = 0
        with opener(ssf, 'rt') as f:
            for l in f:
                p = l.strip().split("\t")
                p[14] = float(p[14])
                read_count[sid] += p[14]
                if not p[9]:
                    continue
                ss_read_count[sid] += p[14]
                if p[9] not in cls:
                    cls[p[9]] = {}
                cls[p[9]][sid] = cls[p[9]].get(sid, 0) + p[14]
                
                if p[10] not in lvl1:
                    lvl1[p[10]] = {}
                lvl1[p[10]][sid] = lvl1[p[10]].get(sid, 0) + p[14]
                
                l1l2 = f"{p[10]}; {p[11]}"
                if l1l2 not in lvl2:
                    lvl2[l1l2] = {}
                lvl2[l1l2][sid] = lvl2[l1l2].get(sid, 0) + p[14]

                if p[12] not in ss:
                    ss[p[12]] = {}
                ss[p[12]][sid] = ss[p[12]].get(sid, 0) + p[14]

                fullss = f"{p[9]}; {p[10]}; {p[11]}; {p[12]}"
                if fullss not in ssall:
                    ssall[fullss] = {}
                ssall[fullss][sid] = ssall[fullss].get(sid, 0) + p[14]



    data = {
        "class" : cls,
        "level1" : lvl1,
        "level2" : lvl2,
        "subsystems" : ss,
        "all_subsystems" : ssall
    }

    os.makedirs(ssdir, exist_ok=True)

    allsids = sorted(read_count.keys())
    
    for what in data:
        ks = sorted(data[what].keys())
        with gzip.open(os.path.join(ssdir, f"{what}_raw.tsv.gz"), "wt") as rawout, gzip.open(os.path.join(ssdir, f"{what}_norm_all_reads.tsv.gz"), "wt") as normallout, gzip.open(os.path.join(ssdir, f"{what}_norm_ss_reads.tsv.gz"), "wt") as normssout:
            print("# Subsystem counts are normalised by dividing by the total number of sequences that matched any protein, and multiplied by 1,000,000", file=normallout)
            print("# Subsystem counts are normalised by dividing by the total number of sequences that matched a protein in a subsystem and multiplied by 1,000,000", file=normssout)
            print("# Subsystem counts are not normalised", file=rawout)
            print("\t".join([what] + allsids), file=rawout)
            print("\t".join([what] + allsids), file=normallout)
            print("\t".join([what] + allsids), file=normssout)
            for ksid in ks:
                for of in rawout, normallout, normssout:
                    print(ksid, file=of, end="")
                for sid in allsids:
                    if sid in data[what][ksid]:
                        print(f"\t{data[what][ksid][sid]}", file=rawout, end="")
                        nc = (data[what][ksid][sid]/read_count[sid]) * 1e6
                        print(f"\t{nc}", file=normallout, end="")
                        ncss = (data[what][ksid][sid]/ss_read_count[sid]) * 1e6
                        print(f"\t{ncss}", file=normssout, end="")
                    else:
                        for of in rawout, normallout, normssout:
                            print("\t0", file=of, end="")
                for of in rawout, normallout, normssout:
                    print(file=of)


try:
    count_ss(snakemake.input.ss, snakemake.output.odir)
except NameError as e:
    print(f"ERROR: {e}", file=sys.stderr)
    parser = argparse.ArgumentParser(description='make file')
    parser.add_argument('-f', help='output file', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()
    make_file(args.f)
except Error as e:
    print(e, file=sys.stderr)

