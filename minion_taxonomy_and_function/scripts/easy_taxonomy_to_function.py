"""
Add the function to the mmseqs easy-taxonomy file.

We want a gzip compressed version of the tophit_report file, and then we look for uniprot IDs in the first column

We use our sqlite database that has the seed to sprot and trembl connections, and the seed subsystems data
"""
import errno
import gzip
import os
import re
import sys
import argparse
import binascii
import sqlite3

__author__ = 'Rob Edwards'

def is_gzip(filename):
    """
    Is this a gzip file?
    """

    """
    This is an elegant solution to test whether a file is gzipped by reading the first two characters.
    I also use a version of this in fastq_pair if you want a C version :)
    See https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed for inspiration
    :param f: the file to test
    :return: True if the file is gzip compressed else false
    """
    with open(filename, 'rb') as i:
        return binascii.hexlify(i.read(2)) == b'1f8b'



def add_ss_tax(tophitfile, sqlitedb, outputfile, verbose=False):
    """
    Add the subsystems to taxonomy
    """

    if not os.path.exists(sqlitedb):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), sqlitedb)

    try:
        con = sqlite3.connect(sqlitedb)
    except sqlite3.Error as e:
        print(f"ERROR Connecting to database: {sqlitedb}\n{e}",
             file=sys.stderr)
        sys.exit(-1)

    urs = re.compile(r'^UniRef\d+_(\w+)')
    cur = con.cursor()

    opener = open
    if is_gzip(tophitfile):
        opener=gzip.open
    
    writeopener = open
    if outputfile.endswith('.gz'):
        writeopener = gzip.open

    #with gzip.open(tophitfile, 'rt') if is_gzip(tophitfile) else open(tophitfile, 'r') as f:
    with opener(tophitfile, 'rt') as f, writeopener(outputfile, 'wt') as out:
        for l in f:
            p = l.strip().split("\t")
            m = urs.match(p[0])
            if m and m.group(1):
                try:
                    # cur.execute("select distinct superclass, class, subclass, subsystem_name, func from subsystems
                    # where func in (select func from trembl where uniprot = ?);", [m.group(1)])
                    cur.execute("select func from trembl where uniprot = ?", [m.group(1)])
                except sqlite3.OperationalError as e:
                    print(f"ERROR: {e} while getting func {m.group(1)}", file=sys.stderr)
                    sys.exit()

                funcs = cur.fetchone()
                if funcs:
                    func = funcs[0]
                    try:
                        cur.execute("select distinct superclass, class, subclass, subsystem_name, func from "
                                    "subsystems where func = ?", [func])
                    except sqlite3.OperationalError as e:
                        print(f"ERROR: {e}\nWhile getting function {func}", file=sys.stderr)
                        sys.exit()
                    counts = 0
                    results = []
                    for s in (cur.fetchall()):
                        results.append("\t".join(list(p) + [func] + list(s)))
                        counts += 1
                    if results:
                        for r in results:
                            print(f"{r}\t{int(p[1])/counts}", file=out)
                    else:
                        print("\t".join(p + [func, "", "", "", "", "", p[1]]), file=out)
                else:
                    print("\t".join(p + ["", "", "", "", "", "", p[1]]), file=out)

            else:
                print(f"Can't parse ID from {p[0]}", file=sys.stderr)
                print("\t".join(p + ["", "", "", "", "", "", p[1]]), file=out)

    con.close()

add_ss_tax(snakemake.input.thr, snakemake.input.sql, snakemake.output.ss, False)
