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


def uniref_to_func(uniref_list):
    """
    Given a list of uniref ids, return the functions of each
    """

    placeholders = ",".join(["?"] * len(uniref_list))  # Create ?, ?, ?
    query = f"SELECT uniprot, func FROM trembl WHERE uniprot IN ({placeholders})"
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


def functions_to_subsystems(functions, known_subsystems, logger=None):
    """
    Given a list of functions, return the subsystems for each
    """

    placeholders = ",".join(["?"] * len(functions))  # Create ?, ?, ?
    query = f"SELECT distinct superclass, class, subclass, subsystem_name, func FROM subsystems WHERE func IN ({placeholders})"
    try:
        cur.execute(query, functions)
    except sqlite3.OperationalError as sqlerr:
        sys.stderr.write("{}".format(sqlerr))
        sys.stderr.write("\nWhile insert on: {}\n".format(p))
        sys.exit()

    for sur in cur.fetchall():
        if sur[4] in known_subsystems:
            known_subsystems[sur[4]].append(list(sur[:4]))
        else:
            known_subsystems[sur[4]] = [list(sur[:4])]
    
    if logger:
        logger.debug(known_subsystems)
    return known_subsystems

def get_taxonomy(list_ids, known_taxonomies, logger):
    """
    Get the taxonomy for a list of ids
    """
    tax_res = pytaxonkit.lineage(list_ids, formatstr="d__{domain|acellular root|superkingdom};p__{phylum};c__{class};o__{order};f__{family};g__{genus};s__{species}")
    for i, t in zip(list_ids, tax_res['Lineage']):
        if not isinstance(t, str):
            t = str(t)
            logger.debug(f"Converted {t} to str for {i}")
        known_taxonomies[i] = t
    if logger:
        logger.debug(known_taxonomies)
    return known_taxonomies


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='mmseqs easy_taxonomy tophit_report.gz', required=True)
    parser.add_argument('-d', help='sqlite3 database of trembl, sprot, seed',
                        default="dummy_database.sqlite")
    parser.add_argument('-o', '--output', help='output file, if not specified, will print to stdout')
    parser.add_argument('-q', '--max_queries', type=int, default=2000,
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
        logging.basicConfig(filename=args.log, encoding='utf-8', level=args.loglevel.upper(),
                            format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    else:
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

    data = {}
    subsystems_cache = {}
    taxonomy_cache = {}
    logger.info(f"Reading data from {args.f}")
    new_taxonomies = set()

    # if the output file exists we back it up and then later we read it and write it to the new file
    backup_file = ""
    if args.output and os.path.exists(args.output):
        backup_file = args.output + ".bak"
        if os.path.exists(backup_file):
            os.remove(backup_file)
        os.rename(args.output, backup_file)

    # set output to stdout if no output file is specified
    if args.output:
        if args.output.endswith('.gz'):
            output_file = gzip.open(args.output, 'wt')
        else:
            output_file = open(args.output, 'w')
    else:
        output_file = sys.stdout

    # read the backup file and write to the output. We need to remember the first column so we can
    # skip it later
    done = set()
    if args.output and os.path.exists(backup_file):
        logger.info(f"Reading backup file {backup_file} to output")
        with gzip.open(backup_file, 'rt') if backup_file.endswith('.gz') else open(backup_file, 'r') as f:
            for ln in f:
                p = ln.strip().split("\t")
                if len(p) < 8:
                    logger.error(f"Line in backup file {backup_file} is malformed: {ln.strip()}")
                    continue
                # remember the first column in done and print everything
                done.add(p[0])
                print(ln.strip(), file=output_file)

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

            if p[0] in done:
                # we have already processed this line
                continue

            m = urs.match(p[0])
            tax = ""
            if m and m.group(1):
                thisid = m.group(1)
                data[thisid] = {}
                data[thisid]['count'] = int(p[1])
                data[thisid]['results'] = p
                data[thisid]['taxid'] = p[5]
                new_taxonomies.add(data[thisid]['taxid'])

                if len(data.keys()) >= args.max_queries:
                    # get all the taxonomies
                    new_taxonomies = new_taxonomies - set(taxonomy_cache.keys())
                    logger.info(f"Getting taxonomies for {len(new_taxonomies)} tax ids")
                    taxonomy_cache = get_taxonomy(new_taxonomies, taxonomy_cache, logger)

                    # get all the functions
                    logger.info(f"Getting functions for {len(data.keys())} uniref ids")
                    uniref_ids = list(data.keys())
                    fns = uniref_to_func(uniref_ids)
                    new_functions = list(set(fns.values())-set(subsystems_cache.keys()))
                    logger.info(f"Getting subsystems for {len(new_functions)} functions")
                    ss = functions_to_subsystems(new_functions, subsystems_cache, logger)
                    logger.info("Printing data")
                    for k in data:
                        if data[k]['taxid'] not in taxonomy_cache:
                            logger.error(f"Taxonomy not found for {k}")
                            taxonomy_cache[data[k]['taxid']]="k__;p__;c__;o__;f__;g__;s__"
                        fn_ss = ["", "", "", "", ""]
                        if k in fns:
                            fn_ss = [fns[k], "", "", "", ""]
                        if k in fns and fns[k] in ss:
                            frac = data[k]['count'] / len(ss[fns[k]])
                            for s in ss[fns[k]]:
                                output = data[k]['results'] + [fns[k]] + s + [str(frac), taxonomy_cache[data[k]['taxid']]]
                                print("\t".join(output), file=output_file)
                        else:
                            output = data[k]['results'] + fn_ss + [str(data[k]['count']), taxonomy_cache[data[k]['taxid']]]
                            print("\t".join(output), file=output_file)
                    data = {}
            else:
                print(f"Can't parse ID from {p[0]}", file=sys.stderr)
                print("\t".join(p + ["", "", "", "", "", "", p[1], ""]), file=output_file)

    con.close()

    # remove the backup file if it exists
    if args.output and os.path.exists(backup_file):
            os.remove(backup_file)
