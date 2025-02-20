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
import logging

__author__ = 'Rob Edwards'


def uniref_to_func(uniref_list, logger):
    """
    Given a list of uniref ids, return the functions of each
    """

    placeholders = ",".join(["?"] * len(uniref_list))  # Create ?, ?, ?
    query = f"SELECT uniprot, func FROM trembl WHERE uniprot IN ({placeholders})"
    logger.debug("Query: {}".format(query))
    try:
        cur.execute(query, uniref_list)
    except sqlite3.OperationalError as sqlerr:
        sys.stderr.write("{}".format(sqlerr))
        sys.stderr.write("\nWhile insert on: {}\n".format(p))
        sys.exit()

    functions = {}
    for fnr in cur.fetchall():
        functions[fnr[0]] = fnr[1]
    return functions


def functions_to_subsystems(functions, known_subsystems, logger):
    """
    Given a list of functions, return the subsystems for each
    """

    placeholders = ",".join(["?"] * len(functions))  # Create ?, ?, ?
    query = f"SELECT distinct superclass, class, subclass, subsystem_name, func FROM subsystems WHERE func IN ({placeholders})"
    logger.debug("Query: {}".format(query))
    try:
        cur.execute(query, functions)
    except sqlite3.OperationalError as sqlerr:
        sys.stderr.write("{}".format(sqlerr))
        sys.stderr.write("\nWhile insert on: {}\n".format(p))
        sys.exit()

    for sur in cur.fetchall():
        if sur[4] in known_subsystems:
            known_subsystems[sur[4]].append(sur[:4])
        else:
            known_subsystems[sur[4]] = [sur[:4]]
    return known_subsystems


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='mmseqs easy_taxonomy tophit_report.gz', required=True)
    parser.add_argument('-d', help='sqlite3 database of trembl, sprot, seed',
                        default="dummy_database.sqlite")
    parser.add_argument('-q', '--max_queries', type=int, default=10000,
                        help='maximum number of simultaneous queries (more requires more memory!)',)
    parser.add_argument('-l', '--log', help='log file')
    parser.add_argument('--loglevel', help='log level', default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'])
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if not os.path.exists(args.d):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.d)

    logger = logging.getLogger(__name__)
    if args.log:
        logging.basicConfig(filename=args.log, encoding='utf-8', level=args.loglevel.upper())
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    try:
        con = sqlite3.connect(args.d)
    except sqlite3.Error as e:
        sys.stderr.write(f"ERROR Connecting to database: {args.d}\n")
        sys.stderr.write(e)
        sys.exit(-1)

    urs = re.compile(r'^UniRef\d+_(\w+)')
    cur = con.cursor()
    logger.info("Connected to database")

    # suggested by chatty, these PRAGMA commands should speed up the database
    cur.execute("PRAGMA cache_size=-8000000;")  # Increase cache size. This is about 8GB
    cur.execute("PRAGMA temp_store=MEMORY;")  # Store temp data in RAM

    taxonomy_cache = {}

    data = {}
    known_subsystems = {}
    logger.info(f"Reading data from {args.f}")
    with gzip.open(args.f, 'rt') if args.f.endswith('.gz') else open(args.f, 'r') as f:
        for ln in f:
            p = ln.strip().split("\t")
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
            tax = ""
            if m and m.group(1):
                thisid = m.group(1)
                data[thisid] = {}
                data[thisid]['count'] = int(p[1])
                data[thisid]['results'] = p

                # get the taxonomy. IN my example test case (above) we often see the same
                # taxon repeated. We cache this, for speed
                logger.debug(f"Getting taxonomy for {p[5]}")
                if p[5] in taxonomy_cache:
                    data[thisid]['tax'] = taxonomy_cache[p[5]]
                else:
                    tax_res = pytaxonkit.lineage([p[5]], prefix=True)
                    tax = tax_res['Lineage'][0]
                    if not isinstance(tax, str):
                        tax = str(tax)
                        print(f"Converted {tax} to str for {p[5]}", file=sys.stderr)
                    taxonomy_cache[p[5]] = tax
                    data[thisid]['tax'] = tax

                if len(data.keys()) % int(args.max_queries/10) == 0:
                    logger.info(f"Processed {len(data.keys())} uniref ids")

                if len(data.keys()) >= args.max_queries:
                    # get all the functions
                    logger.info(f"Getting functions for {len(data.keys())} uniref ids")
                    uniref_ids = list(data.keys())
                    fns = uniref_to_func(uniref_ids, logger)
                    new_functions = set(fns.values())-set(known_subsystems.keys())
                    ss = functions_to_subsystems(new_functions, known_subsystems, logger)
                    for k in data:
                        fn_ss = ["", "", "", "", ""]
                        if k in fns:
                            fn_ss = [fns[k], "", "", "", ""]
                        if fns[k] in ss:
                            frac = data[k]['count'] / len(ss[fns[k]])
                            for s in ss[fns[k]]:
                                print(f"{data[k]['results']}\t{fns[k]}\t{s}\t{frac}\t{data[k]['tax']}")
                        else:
                            print(f"{data[k]['results']}\t{fn_ss}\t{data[k]['count']}\t{data[k]['tax']}")
                    data = {}
                    num_lines_processed = 0
            else:
                print(f"Can't parse ID from {p[0]}", file=sys.stderr)
                print("\t".join(p + ["", "", "", "", "", "", p[1], ""]))

    con.close()
