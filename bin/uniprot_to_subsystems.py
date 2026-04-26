"""
Download the latest BV-BRC subsystems and connect them to UniProt IDs.

This script:
  1. Downloads the patric_uniprot_linkout.gz mapping file from the BV-BRC FTP
     server.  The file has two tab-separated columns: UniProt ID and patric_id
     (e.g. "fig|511145.12.peg.380").
  2. Extracts every unique genome ID that is encoded in the patric_id values.
  3. For each genome ID, checks whether a PATRIC subsystem annotation file
     exists on the FTP server and, if so, downloads it.  Downloads are run in
     parallel (using multiple processes) with a configurable number of workers
     and a per-request delay so as not to overload the server.
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

Checkpoint / resume
-------------------
After each genome is downloaded the result is appended to an intermediate file
(--intermediate, default: subsystems_intermediate.tsv).  A companion
"<intermediate>.done" file records every genome ID that has been attempted
(including genomes that had no subsystem file).  On restart the script reads
both files, skips genomes that are already done, and merges the cached data
with any new downloads, so a run that was interrupted can be completed without
re-downloading data that was already fetched.

Multiple instances
------------------
Genome IDs are shuffled into a random order before downloading (use
--no-randomize to disable).  Running several instances simultaneously will
therefore tend to download disjoint sets of genomes, reducing redundant work.
Be aware that running many instances at once increases the risk of overloading
the BV-BRC FTP server; always use a reasonable --delay value.
"""

import ftplib
import gzip
import io
import os
import random
import sys
import argparse
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

__author__ = 'Rob Edwards'

FTP_HOST = 'ftp.bv-brc.org'
LINKOUT_PATH = '/linkouts/uniprot/patric_uniprot_linkout.gz'
SUBSYSTEM_PATH_TEMPLATE = '/genomes/{genome_id}/{genome_id}.PATRIC.subsystem.tab'


_FTP_RETRIES = 3          # number of attempts before giving up
_FTP_RETRY_DELAY = 5.0   # seconds to wait between retry attempts


class ReusableFTP_TLS(ftplib.FTP_TLS):
    """
        Explicit FTPS with shared TLS session support.

        This is required, because we are using ftps
    """
    def ntransfercmd(self, cmd, rest=None):
        conn, size = super().ntransfercmd(cmd, rest)
        if self._prot_p:
            # Wrap the data connection socket with the existing control channel session
            conn = self.context.wrap_socket(conn, 
                                            server_hostname=self.host, 
                                            session=self.sock.session)
        return conn, size


def _connect():
    """Open and return an anonymous FTP connection to FTP_HOST.

    Note: We use FTP_TLS for the ftps:// connection.
    """
    ftps = ReusableFTP_TLS(FTP_HOST, timeout=60)
    ftps.login('anonymous', 'guest')
    ftps.prot_p()
    return ftps


def _ftp_retr(path):
    """Download *path* from the BV-BRC FTP server with retry logic.

    Parameters
    ----------
    path : str
        Absolute path on the FTP server.

    Returns
    -------
    io.BytesIO
        Buffer containing the raw file bytes.

    Raises
    ------
    ftplib.error_perm
        Re-raised if the server returns a permanent error (e.g. 550 file not
        found) — callers interpret this as "file does not exist".
    Exception
        Re-raised after *_FTP_RETRIES* failed attempts for any other error.
    """
    last_exc = None
    for attempt in range(1, _FTP_RETRIES + 1):
        try:
            ftp = _connect()
            buf = io.BytesIO()
            ftp.retrbinary(f'RETR {path}', buf.write)
            ftp.quit()
            buf.seek(0)
            return buf
        except ftplib.error_perm:
            # Permanent error (e.g. 550 – file not found); do not retry.
            raise
        except Exception as exc:
            last_exc = exc
            if attempt < _FTP_RETRIES:
                time.sleep(_FTP_RETRY_DELAY)
    raise last_exc


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
        try:
            buf = _ftp_retr(LINKOUT_PATH)
        except Exception as exc:
            print(
                f"Error: could not download linkout file from {FTP_HOST}{LINKOUT_PATH}: {exc}\n"
                f"Hint: if you have already downloaded the file, pass it with -l / --linkout "
                f"to skip the FTP step.",
                file=sys.stderr,
            )
            sys.exit(1)
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
        buf = _ftp_retr(path)
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


def load_intermediate(intermediate_path, verbose=False):
    """Load already-completed genomes from checkpoint files.

    Two companion files are read:

    * ``<intermediate_path>.done`` – one genome ID per line; lists every genome
      that has already been attempted (including genomes that had no subsystem
      file, so they are not re-downloaded on restart).
    * ``<intermediate_path>`` – tab-separated data written by previous runs:
      patric_id, role_name, superclass, class, subclass, subsystem_name.

    Parameters
    ----------
    intermediate_path : str
    verbose : bool

    Returns
    -------
    done_genomes : set of str
    cached_data : dict
        patric_id -> (role_name, superclass, class, subclass, subsystem_name)
    """
    done_genomes = set()
    cached_data = {}
    done_path = intermediate_path + '.done'

    if os.path.exists(done_path):
        with open(done_path, 'r') as fh:
            for line in fh:
                line = line.strip()
                if line:
                    done_genomes.add(line)
        if verbose:
            print(
                f"Loaded {len(done_genomes):,} already-processed genome IDs "
                f"from {done_path}",
                file=sys.stderr,
            )

    if os.path.exists(intermediate_path):
        open_func = gzip.open if intermediate_path.endswith('.gz') else open
        with open_func(intermediate_path, 'rt') as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 6:
                    continue
                # format: patric_id, role_name, superclass, class, subclass,
                #         subsystem_name
                patric_id = parts[0]
                cached_data[patric_id] = (
                    parts[1], parts[2], parts[3], parts[4], parts[5]
                )
        if verbose:
            print(
                f"Loaded {len(cached_data):,} cached patric entries "
                f"from {intermediate_path}",
                file=sys.stderr,
            )

    return done_genomes, cached_data


def _append_to_intermediate(intermediate_path, done_path, genome_id, result):
    """Append a single genome's results to the checkpoint files.

    Parameters
    ----------
    intermediate_path : str
    done_path : str
    genome_id : str
    result : dict or None
        patric_id -> (role_name, superclass, class, subclass, subsystem_name)
    """
    if result:
        open_func = gzip.open if intermediate_path.endswith('.gz') else open
        with open_func(intermediate_path, 'at') as fh:
            for patric_id, fields in result.items():
                fh.write(patric_id + '\t' + '\t'.join(fields) + '\n')
    with open(done_path, 'a') as fh:
        fh.write(genome_id + '\n')


def download_subsystems(
    genome_ids,
    workers=5,
    delay=0.5,
    intermediate_path=None,
    randomize=True,
    verbose=False,
):
    """Download subsystem files for all genome IDs in parallel.

    Uses ``ProcessPoolExecutor`` to exploit multiple CPU cores on HPC systems.
    Each worker process manages its own FTP connection, so connections are not
    shared across processes.

    After each genome is downloaded the result is written to *intermediate_path*
    so that a run that is interrupted can be restarted without re-downloading
    data that was already fetched.

    Parameters
    ----------
    genome_ids : iterable of str
    workers : int
        Maximum number of parallel worker processes.
    delay : float
        Per-request sleep (seconds) inside each worker to rate-limit requests.
    intermediate_path : str or None
        Path to the checkpoint data file.  Pass ``None`` to disable checkpointing.
    randomize : bool
        Shuffle the genome list before processing.  Useful when running
        multiple simultaneous instances so they tend to work on different
        genomes.  Note that running many concurrent instances with a short
        *delay* increases the risk of overloading the BV-BRC FTP server.
    verbose : bool

    Returns
    -------
    subsystem_data : dict
        Mapping patric_id -> (role_name, superclass, class, subclass,
        subsystem_name) for all genomes that had a subsystem file.
    """
    # --- load checkpoint ---
    done_genomes = set()
    subsystem_data = {}
    done_path = intermediate_path + '.done' if intermediate_path else None

    if intermediate_path:
        done_genomes, subsystem_data = load_intermediate(
            intermediate_path, verbose=verbose
        )

    genome_list = [g for g in genome_ids if g not in done_genomes]

    if randomize:
        random.shuffle(genome_list)

    total = len(genome_list)
    skipped = len(done_genomes)

    if verbose:
        print(
            f"Downloading subsystem files for {total:,} genomes "
            f"({skipped:,} already done) using {workers} worker(s) "
            f"with {delay}s delay...",
            file=sys.stderr,
        )

    completed = 0
    found = 0

    with ProcessPoolExecutor(max_workers=workers) as executor:
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

            if intermediate_path:
                _append_to_intermediate(
                    intermediate_path, done_path, genome_id, result
                )

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
        '-t', '--workers',
        help='number of parallel worker processes [Default: %(default)s]',
        type=int,
        default=5,
    )
    parser.add_argument(
        '-d', '--delay',
        help=(
            'delay in seconds between FTP requests per worker to avoid '
            'overloading the server [Default: %(default)s]'
        ),
        type=float,
        default=0.5,
    )
    parser.add_argument(
        '-i', '--intermediate',
        help=(
            'intermediate checkpoint file; already-downloaded data is appended '
            'here so a crashed run can be resumed. A companion '
            '"<file>.done" file tracks completed genome IDs. '
            '[Default: %(default)s]'
        ),
        default='subsystems_intermediate.tsv',
    )
    parser.add_argument(
        '--no-randomize',
        help=(
            'process genomes in sorted order instead of shuffling. '
            'Shuffling (the default) allows multiple simultaneous instances '
            'to work on different genomes; disable when reproducible ordering '
            'is required.'
        ),
        action='store_true',
        default=False,
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
        workers=args.workers,
        delay=args.delay,
        intermediate_path=args.intermediate,
        randomize=not args.no_randomize,
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
