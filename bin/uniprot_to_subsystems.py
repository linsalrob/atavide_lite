"""
Download the latest BV-BRC subsystems and connect them to UniProt IDs.

This script:
  1. Downloads the patric_uniprot_linkout.gz mapping file from the BV-BRC FTP
     server.  The file has two tab-separated columns: UniProt ID and patric_id
     (e.g. "fig|511145.12.peg.380").
  2. Extracts every unique genome ID that is encoded in the patric_id values.
  3. For each genome ID, checks whether a PATRIC subsystem annotation file
     exists on the FTP server and, if so, downloads it.  Downloads are run in
     parallel with a configurable number of threads and a per-request delay so
     as not to overload the server.
  4. Merges the two data sources on patric_id and writes a tab-separated output
     file with the following columns:

        0  uniprot_id
        1  patric_id
        2  role_name
        3  superclass
        4  class
        5  subclass
        6  subsystem_name

Output can be written as plain text or gzip-compressed (use a .gz extension).
"""

import ftplib
import gzip
import io
import os
import sys
import argparse
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

__author__ = 'Rob Edwards'

FTP_HOST = 'ftp.bv-brc.org'
LINKOUT_PATH = '/linkouts/uniprot/patric_uniprot_linkout.gz'
SUBSYSTEM_PATH_TEMPLATE = '/genomes/{genome_id}/{genome_id}.PATRIC.subsystem.tab'


def _connect():
    """Open and return an anonymous FTP connection to FTP_HOST."""
    ftp = ftplib.FTP(FTP_HOST, timeout=60)
    ftp.login()
    return ftp


def download_linkout(linkout_file=None, verbose=False):
    """Download (or read from disk) and parse the patric_uniprot_linkout file.

    Parameters
    ----------
    linkout_file : str or None
        If given, read from this local file instead of downloading.
    verbose : bool

    Returns
    -------
    patric_to_uniprot : dict
        Mapping patric_id -> uniprot_id.
    genome_ids : set
        Set of genome IDs extracted from the patric_ids.
    """
    if linkout_file:
        if verbose:
            print(f"Reading linkout file from {linkout_file}", file=sys.stderr)
        opener = gzip.open(linkout_file, 'rt') if linkout_file.endswith('.gz') else open(linkout_file, 'r')
    else:
        if verbose:
            print(f"Connecting to {FTP_HOST} to download linkout file...", file=sys.stderr)
        ftp = _connect()
        buf = io.BytesIO()
        ftp.retrbinary(f'RETR {LINKOUT_PATH}', buf.write)
        ftp.quit()
        buf.seek(0)
        opener = gzip.open(buf, 'rt')

    patric_to_uniprot = {}
    genome_ids = set()

    with opener as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            uniprot_id = parts[0]
            patric_id = parts[1]
            # patric_id format: fig|{genome_id}.peg.{number}
            if patric_id.startswith('fig|'):
                inner = patric_id[4:]   # strip the 'fig|' prefix
                dot_peg = inner.find('.peg.')
                if dot_peg != -1:
                    genome_ids.add(inner[:dot_peg])
            patric_to_uniprot[patric_id] = uniprot_id

    if verbose:
        print(
            f"Found {len(patric_to_uniprot):,} patric->uniprot mappings "
            f"across {len(genome_ids):,} unique genomes",
            file=sys.stderr,
        )

    return patric_to_uniprot, genome_ids


def _download_one_genome(genome_id, delay):
    """Download and parse the subsystem file for a single genome.

    Parameters
    ----------
    genome_id : str
    delay : float
        Seconds to sleep before making the FTP request (rate limiting).

    Returns
    -------
    (genome_id, result)
        result is a dict mapping patric_id ->
        (role_name, superclass, class, subclass, subsystem_name),
        or None if the file does not exist or an error occurred.
    """
    if delay > 0:
        time.sleep(delay)

    path = SUBSYSTEM_PATH_TEMPLATE.format(genome_id=genome_id)

    try:
        ftp = _connect()
        buf = io.BytesIO()
        ftp.retrbinary(f'RETR {path}', buf.write)
        ftp.quit()
    except ftplib.error_perm:
        # 550 – file does not exist; this is expected for many genomes
        return genome_id, None
    except Exception as exc:
        sys.stderr.write(f"Warning: could not download {path}: {exc}\n")
        return genome_id, None

    result = {}
    for line in buf.getvalue().decode('utf-8', errors='replace').splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) < 11:
            continue
        patric_id    = parts[2]
        role_name    = parts[6]
        superclass   = parts[7]
        cls          = parts[8]
        subclass     = parts[9]
        subsystem_nm = parts[10]
        result[patric_id] = (role_name, superclass, cls, subclass, subsystem_nm)

    return genome_id, result


def download_subsystems(genome_ids, threads=5, delay=0.5, verbose=False):
    """Download subsystem files for all genome IDs in parallel.

    Parameters
    ----------
    genome_ids : iterable of str
    threads : int
        Maximum number of concurrent FTP connections.
    delay : float
        Per-thread delay (seconds) between requests.
    verbose : bool

    Returns
    -------
    subsystem_data : dict
        Mapping patric_id -> (role_name, superclass, class, subclass,
        subsystem_name) for all genomes that had a subsystem file.
    """
    genome_list = sorted(genome_ids)
    total = len(genome_list)

    if verbose:
        print(
            f"Downloading subsystem files for {total:,} genomes "
            f"using {threads} thread(s) with {delay}s delay...",
            file=sys.stderr,
        )

    subsystem_data = {}
    completed = 0
    found = 0

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(_download_one_genome, gid, delay): gid
            for gid in genome_list
        }
        for future in as_completed(futures):
            genome_id, result = future.result()
            completed += 1
            if result is not None:
                subsystem_data.update(result)
                found += 1
                if verbose and found % 100 == 0:
                    print(
                        f"  [{completed}/{total}] Downloaded subsystems for "
                        f"{found} genomes so far...",
                        file=sys.stderr,
                    )
            elif verbose and completed % 500 == 0:
                print(f"  [{completed}/{total}] processed...", file=sys.stderr)

    if verbose:
        print(
            f"Finished: {found:,} genomes had subsystem files; "
            f"{len(subsystem_data):,} patric entries collected",
            file=sys.stderr,
        )

    return subsystem_data


def write_output(patric_to_uniprot, subsystem_data, output_path, verbose=False):
    """Merge the two datasets and write the output file.

    Parameters
    ----------
    patric_to_uniprot : dict
        patric_id -> uniprot_id
    subsystem_data : dict
        patric_id -> (role_name, superclass, class, subclass, subsystem_name)
    output_path : str
    verbose : bool
    """
    open_func = gzip.open if output_path.endswith('.gz') else open

    written = 0
    with open_func(output_path, 'wt') as out:
        out.write("uniprot_id\tpatric_id\trole_name\tsuperclass\tclass\tsubclass\tsubsystem_name\n")
        for patric_id in sorted(patric_to_uniprot):
            if patric_id not in subsystem_data:
                continue
            uniprot_id = patric_to_uniprot[patric_id]
            role_name, superclass, cls, subclass, subsystem_nm = subsystem_data[patric_id]
            out.write(
                f"{uniprot_id}\t{patric_id}\t{role_name}\t"
                f"{superclass}\t{cls}\t{subclass}\t{subsystem_nm}\n"
            )
            written += 1

    if verbose:
        print(f"Wrote {written:,} merged entries to {output_path}", file=sys.stderr)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=(
            'Download BV-BRC subsystem annotations and map them to UniProt IDs. '
            'Produces a tab-separated file with columns: '
            'uniprot_id, patric_id, role_name, superclass, class, subclass, subsystem_name.'
        )
    )
    parser.add_argument(
        '-o', '--output',
        help='output file; use a .gz extension for gzip compression [Default: %(default)s]',
        default='uniprot_to_subsystems.tsv.gz',
    )
    parser.add_argument(
        '-l', '--linkout',
        help=(
            'use a pre-downloaded patric_uniprot_linkout.gz file instead of '
            'fetching it from the FTP server'
        ),
        default=None,
    )
    parser.add_argument(
        '-t', '--threads',
        help='number of parallel download threads [Default: %(default)s]',
        type=int,
        default=5,
    )
    parser.add_argument(
        '-d', '--delay',
        help=(
            'delay in seconds between FTP requests per thread to avoid '
            'overloading the server [Default: %(default)s]'
        ),
        type=float,
        default=0.5,
    )
    parser.add_argument(
        '-v', '--verbose',
        help='verbose output',
        action='store_true',
    )
    args = parser.parse_args()

    # Step 1: obtain the UniProt <-> patric_id mapping
    patric_to_uniprot, genome_ids = download_linkout(
        linkout_file=args.linkout,
        verbose=args.verbose,
    )

    if not genome_ids:
        print("No genome IDs found in the linkout file. Exiting.", file=sys.stderr)
        sys.exit(1)

    # Step 2: download subsystem annotation files in parallel
    subsystem_data = download_subsystems(
        genome_ids,
        threads=args.threads,
        delay=args.delay,
        verbose=args.verbose,
    )

    if not subsystem_data:
        print("No subsystem data could be downloaded. Exiting.", file=sys.stderr)
        sys.exit(1)

    # Step 3: merge and write output
    write_output(
        patric_to_uniprot,
        subsystem_data,
        args.output,
        verbose=args.verbose,
    )

    if args.verbose:
        print("Done.", file=sys.stderr)
