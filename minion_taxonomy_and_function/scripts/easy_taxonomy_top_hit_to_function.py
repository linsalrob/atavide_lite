"""
`mmseqs easy-taxonomy` reports two key files:
1. `_tophit_aln.gz` has the top alignment(s) for the sequence. If there are equally valid alignments, both will be reported
2. `_tophit_report.gz` has the protein that matches and the number of reads that match to it.

The `easy_taxonomy_to_function' adds the subsystems to the first second file. This will add the subsystems to the
first version, with a few interactive choices!
"""

import os
import random
import sys
import argparse
import errno
import gzip
import re
import binascii
import sqlite3

__author__ = 'Rob Edwards'


def is_gzip(filename):
    """
    Is this a gzip file?
    """
    with open(filename, 'rb') as i:
        return binascii.hexlify(i.read(2)) == b'1f8b'

def printout(results, mthd, outputfile):
    """
    print the data depending on the method:
    'f': print the first subsystem result
    'r': print one at random
    'a': print all results
    's': print all subsystems
    't': the first hit
    """

    if mthd == 'r':
        print(random.choice(results['subsystems'] + results['nosubsystems']), file=outputfile)
    if mthd == 's':
        for r in results['subsystems']:
            print(r, file=outputfile)
    if mthd == 'a':
        for r in results['subsystems'] + results['nosubsystems']:
            print(r, file=outputfile)
    if mthd == 'f':
        if results['subsystems']:
            print(results['subsystems'][0], file=outputfile)
    if mthd == 't':
        if results['tophit']:
            print(results['tophit'], file=outputfile)

def add_ss_tax(tophitaln, mthd, sqlitedb, outputfile, verbose=False):
    """
    Add the subsystems to the top hit alignment. We use the method in
    mthd to choose whether to keep all alignments, a random alignment,
    the one with the subsystem
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

    opener = gzip.open if is_gzip(tophitaln) else open
    writeopener = gzip.open if outputfile.endswith('.gz') else open

    with opener(tophitaln, 'rt') as f, writeopener(outputfile, 'wt') as out:
        lastsequence = None
        results = {'subsystems': [], 'nosubsystems': [], 'tophit': None}
        for l in f:
            p = l.strip().split("\t")
            if not lastsequence:
                lastsequence = p[0]
            if lastsequence != p[0]:
                printout(results, mthd, out)
                lastsequence = p[0]
                results = {'subsystems': [], 'nosubsystems': [], 'tophit': None}
            m = urs.match(p[1])
            if m and m.group(1):
                try:
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
                    for s in (cur.fetchall()):
                        results['subsystems'].append("\t".join(list(p) + [func] + list(s)))
                        if not results['tophit']:
                            results['tophit'] = "\t".join(list(p) + [func] + list(s))
                        counts += 1
                    if counts == 0:
                        results['nosubsystems'].append("\t".join(list(p) + [func] + ["", "", "", "", ""]))
                        if not results['tophit']:
                            results['tophit'] = "\t".join(list(p) + [func] + ["", "", "", "", ""])
                else:
                    results['nosubsystems'].append("\t".join(list(p) + ["", "", "", "", "", ""]))
                    if not results['tophit']:
                        results['tophit'] = "\t".join(list(p) + ["", "", "", "", "", ""])
            else:
                print(f"Can't parse ID from {p[0]}", file=sys.stderr)
                results['nosubsystems'].append("\t".join(list(p) + ["", "", "", "", "", ""]))
                if not results['tophit']:
                    results['tophit'] = "\t".join(list(p) + ["", "", "", "", "", ""])



    con.close()


if __name__ == "__main__":
    valid_methods = ['f', 'a', 'r', 's', 't']
    parser = argparse.ArgumentParser(description='Print some ')
    parser.add_argument('-a', '--tophitaln',  help='mmseqs top hit alignment file', required=True)
    parser.add_argument('-o', '--outputfile', help='input file', required=True)
    parser.add_argument('-m', '--method', choices=valid_methods, required=True,
                        help='method to choose the print out. Currrent choices are\n' +
                        'f: print out the first subsystem match\n' +
                        'r: print one match at random from subsystems + no subsytems\n' +
                        'a: print all matches\ns: print all subsystem matches\n' +
                        't: use the top hit in the alignment file')
    parser.add_argument('-s', '--sqlitedb', help='sqlite database', required=True)
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    add_ss_tax(args.tophitaln, args.method, args.sqlitedb, args.outputfile, args.verbose)
