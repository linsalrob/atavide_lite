"""
Add the function and the full taxa to the mmseqs easy-taxonomy file.

We want a gzip compressed version of the tophit_report file, and then we look for uniprot IDs in the first column

We use our sqlite database that has the seed to sprot and trembl connections, and the seed subsystems data. We also use
"""
import errno
import gzip
import os
import re
import sys
import argparse
import sqlite3
import pytaxonkit

__author__ = 'Rob Edwards'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='mmseqs easy_taxonomy tophit_report.gz', required=True)
    parser.add_argument('-d', help='sqlite3 database of trembl, sprot, seed', default="dummy_database.sqlite")
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if not os.path.exists(args.d):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.d)

    try:
        con = sqlite3.connect(args.d)
    except sqlite3.Error as e:
        sys.stderr.write(f"ERROR Connecting to database: {args.d}\n")
        sys.stderr.write(e)
        sys.exit(-1)

    urs = re.compile(r'^UniRef\d+_(\w+)')
    cur = con.cursor()

    taxonomy_cache = {}

    with gzip.open(args.f, 'rt') if args.f.endswith('.gz') else open(args.f, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            # columns are
            # (0) Target identifier uniref id eg: UniRef50_S9RHF5
            # (1) Number of sequences aligning to target eg: 21
            # (2) Unique coverage of target uniqueAlignedResidues / totalLength eg: 0.988
            # (3) Target coverage alignedResidues / totalLength: eg 20.756
            # (4) Average sequence identity: eg 0.433
            # (5) taxonomy id eg: 1077464
            # (6) rank eg: subspecies
            # (7) scientific name: Streptococcus oralis subsp.tigurinus
            m = urs.match(p[0])
            if m and m.group(1):

                # get the taxonomy. IN my example test case (above) we often see the same
                # taxon repeated. We cache this, for speed
                tax = ""
                if p[5] in taxonomy_cache:
                    tax = taxonomy_cache[p[5]]
                else:
                    tax_res = pytaxonkit.lineage(p[5])
                    taxonomy_cache[p[5]] = tax_res['Lineage'][0]
                    tax = tax_res['Lineage'][0]

                try:
                    # cur.execute("select distinct superclass, class, subclass, subsystem_name, func from subsystems
                    # where func in (select func from trembl where uniprot = ?);", [m.group(1)])
                    cur.execute("select func from trembl where uniprot = ?", [m.group(1)])
                except sqlite3.OperationalError as e:
                    sys.stderr.write("{}".format(e))
                    sys.stderr.write("\nWhile insert on: {}\n".format(p))
                    sys.exit()

                funcs = cur.fetchone()
                if funcs:
                    func = funcs[0]
                    try:
                        cur.execute("select distinct superclass, class, subclass, subsystem_name, func from "
                                    "subsystems where func = ?", [func])
                    except sqlite3.OperationalError as e:
                        sys.stderr.write("{}".format(e))
                        sys.stderr.write("\nWhile insert on: {}\n".format(p))
                        sys.exit()
                    counts = 0
                    results = []
                    for s in (cur.fetchall()):
                        # print(">>", end="", file=sys.stderr)
                        # print("|\t|".join(s), end="", file=sys.stderr);
                        # print("<<", file=sys.stderr)
                        results.append("\t".join(list(p) + [func] + list(s)))
                        counts += 1
                    if results:
                        for r in results:
                            print(f"{r}\t{int(p[1])/counts}\t{tax}")
                    else:
                        # print(f"Can't find a class for {m.group(1)}", file=sys.stderr)
                        print("\t".join(p + [func, "", "", "", "", "", p[1], tax]))
                else:
                    print("\t".join(p + ["", "", "", "", "", "", p[1], tax]))

            else:
                print(f"Can't parse ID from {p[0]}", file=sys.stderr)
                print("\t".join(p + ["", "", "", "", "", "", p[1], ""]))

    con.close()
