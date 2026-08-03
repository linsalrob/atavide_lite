#!/usr/bin/env python3
"""Build a normalized SQLite database mapping UniProt IDs to BV-BRC subsystem annotations.

Input files
-----------
1. patric_uniprot_linkout.gz (or plain .tsv)
     Tab-separated, no header:  uniprot_id <TAB> patric_id
     Produced by: uniprot_to_subsystems.py  (-s / --save-linkout)

2. Subsystem intermediate TSV  (plain or .gz)
     Tab-separated, no header:
       patric_id, role_name, superclass, class, subclass, subsystem_name
     Produced by: uniprot_to_subsystems.py  (-i / --intermediate)

Normalized schema
-----------------
  roles              (role_id, role_name)
  subsystems         (subsystem_id, superclass, class, subclass, subsystem_name)
  role_subsystems    (role_id, subsystem_id)
  uniprot_patric_role(uniprot_id, patric_id, role_id)

Key relationships
-----------------
* Many patric_ids can map to the same role_name.
* A role_name can appear in more than one subsystem.
* subsystem_name is globally unique.
* Not every UniProt/PATRIC mapping will necessarily map to a subsystem.

Staging tables (dropped unless --keep-staging)
----------------------------------------------
  patric_role_staging(patric_id, role_id)  -- patric_id → role_id
  linkout_staging    (uniprot_id, patric_id) -- raw linkout rows

Performance notes
-----------------
* SQLite import pragmas are applied before any data load.
* All indexes are created AFTER the data is loaded (not before).
* Staging indexes are created after staging tables are populated, and before
  the final JOIN that populates uniprot_patric_role.
* Data is inserted in batches (--batch-size, default 50 000) within explicit
  transactions so that a single COMMIT covers many rows.
* Roles and subsystems are small enough to cache in a Python dict; the
  potentially large patric_id→role_id mapping is kept in a SQLite staging
  table to avoid excessive RAM use.

Usage
-----
  python subsystems_to_sqlite.py \\
      --linkout    patric_uniprot_linkout.gz \\
      --subsystems subsystems_intermediate.tsv \\
      --db         subsystems.db

Dependencies
------------
Python 3 standard library only: sqlite3, csv, argparse, gzip, pathlib, sys.
"""

import argparse
import csv
import gzip
import pathlib
import sqlite3
import sys

__author__ = 'Rob Edwards'

DEFAULT_BATCH_SIZE = 50_000
PROGRESS_INTERVAL  = 1_000_000   # print a status line every N rows

# ---------------------------------------------------------------------------
# SQLite performance pragmas – applied once at connection open.
# WAL mode allows concurrent reads during the load; NORMAL sync is safe
# for bulk loads (no fsync after every write); MEMORY temp store and a
# large page cache dramatically reduce I/O.
# ---------------------------------------------------------------------------
_PRAGMAS = [
    "PRAGMA journal_mode = WAL",
    "PRAGMA synchronous = NORMAL",
    "PRAGMA temp_store = MEMORY",
    "PRAGMA cache_size = -200000",   # ~200 MB page cache
]

# ---------------------------------------------------------------------------
# Schema: permanent tables
# ---------------------------------------------------------------------------
_DDL_PERMANENT = """\
CREATE TABLE IF NOT EXISTS roles (
    role_id   INTEGER PRIMARY KEY,
    role_name TEXT    NOT NULL UNIQUE
);

CREATE TABLE IF NOT EXISTS subsystems (
    subsystem_id   INTEGER PRIMARY KEY,
    superclass     TEXT,
    class          TEXT,
    subclass       TEXT,
    subsystem_name TEXT NOT NULL UNIQUE
);

-- Many-to-many: one role can appear in multiple subsystems and one
-- subsystem contains multiple roles.
CREATE TABLE IF NOT EXISTS role_subsystems (
    role_id      INTEGER NOT NULL,
    subsystem_id INTEGER NOT NULL,
    PRIMARY KEY (role_id, subsystem_id)
);

-- Main mapping: UniProt protein ID → PATRIC protein ID → functional role.
-- Not every UniProt/PATRIC mapping will have a subsystem annotation.
CREATE TABLE IF NOT EXISTS uniprot_patric_role (
    uniprot_id TEXT    NOT NULL,
    patric_id  TEXT    NOT NULL,
    role_id    INTEGER NOT NULL,
    PRIMARY KEY (uniprot_id, patric_id, role_id)
);
"""

# ---------------------------------------------------------------------------
# Schema: staging tables (populated during load, dropped afterwards unless
# --keep-staging is set).
# We drop these before (re-)creating them so that a restarted run always
# begins with clean staging data.
# ---------------------------------------------------------------------------
_DROP_STAGING = """\
DROP TABLE IF EXISTS patric_role_staging;
DROP TABLE IF EXISTS linkout_staging;
"""

_DDL_STAGING = """\
-- patric_id → role_id mapping derived from the subsystem annotation files.
-- May contain tens of millions of rows; indexed before the final JOIN.
CREATE TABLE patric_role_staging (
    patric_id TEXT    NOT NULL,
    role_id   INTEGER NOT NULL
);

-- Raw rows from the patric_uniprot_linkout file.
CREATE TABLE linkout_staging (
    uniprot_id TEXT NOT NULL,
    patric_id  TEXT NOT NULL
);
"""

# Staging indexes – created AFTER staging tables are fully populated so that
# bulk inserts are not slowed down by index maintenance.
_INDEX_STAGING = [
    "CREATE INDEX idx_stg_patric_role    ON patric_role_staging(patric_id)",
    "CREATE INDEX idx_stg_linkout_patric ON linkout_staging(patric_id)",
]

# Final indexes on permanent tables – created AFTER uniprot_patric_role is
# populated, then ANALYZE is run so the query planner uses them.
_INDEX_FINAL = [
    "CREATE INDEX IF NOT EXISTS idx_upr_uniprot    ON uniprot_patric_role(uniprot_id)",
    "CREATE INDEX IF NOT EXISTS idx_upr_patric     ON uniprot_patric_role(patric_id)",
    "CREATE INDEX IF NOT EXISTS idx_upr_role       ON uniprot_patric_role(role_id)",
    "CREATE INDEX IF NOT EXISTS idx_rs_subsystem   ON role_subsystems(subsystem_id)",
    "CREATE INDEX IF NOT EXISTS idx_sub_name       ON subsystems(subsystem_name)",
    "CREATE INDEX IF NOT EXISTS idx_sub_hierarchy  ON subsystems(superclass, class, subclass)",
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _open_file(path):
    """Open a plain-text or gzip-compressed file for reading (UTF-8).

    Exits with a clear error message if the file cannot be opened.
    """
    p = pathlib.Path(path)
    if not p.exists():
        print(f"Error: file not found: {path}", file=sys.stderr)
        sys.exit(1)
    try:
        if p.suffix == '.gz':
            return gzip.open(p, 'rt', encoding='utf-8', errors='replace')
        return open(p, 'r', encoding='utf-8', errors='replace')
    except OSError as exc:
        print(f"Error: cannot open {path}: {exc}", file=sys.stderr)
        sys.exit(1)


def _setup_db(db_path):
    """Create (or open) the SQLite database and apply schema and pragmas.

    Permanent tables are created with IF NOT EXISTS so that an existing
    database can be extended.  Staging tables are always dropped and
    recreated to guarantee a clean slate for each run.
    """
    conn = sqlite3.connect(db_path)
    # Disable Python's implicit transaction management; we issue explicit
    # BEGIN / COMMIT for maximum bulk-load performance.
    conn.isolation_level = None

    for pragma in _PRAGMAS:
        conn.execute(pragma)

    # executescript issues an implicit COMMIT before running, which is fine
    # here because no transaction is open yet.
    conn.executescript(_DDL_PERMANENT)
    conn.executescript(_DROP_STAGING)
    conn.executescript(_DDL_STAGING)

    return conn


# ---------------------------------------------------------------------------
# Step 1 – load subsystem annotations
# ---------------------------------------------------------------------------

def load_subsystems(conn, subsystem_path, batch_size, verbose):
    """Stream the subsystem intermediate file and populate:

    * ``roles``
    * ``subsystems``
    * ``role_subsystems``
    * ``patric_role_staging``

    Expected columns (TSV, no header):

    ======  ==============
    col 0   patric_id
    col 1   role_name
    col 2   superclass     (may be empty → stored as NULL)
    col 3   class          (may be empty → stored as NULL)
    col 4   subclass       (may be empty → stored as NULL)
    col 5   subsystem_name (may be empty → row's subsystem link skipped)
    ======  ==============
    """
    # In-memory lookup caches – roles and subsystems number in the thousands,
    # so keeping them in RAM is safe and avoids repeated SELECT queries.
    role_cache      = {}   # role_name      → role_id
    subsystem_cache = {}   # subsystem_name → subsystem_id
    # Track (role_id, subsystem_id) pairs already added to role_subsystems
    # to avoid issuing redundant INSERT OR IGNORE calls within a batch.
    rs_seen = set()

    rs_batch      = []   # pending role_subsystems rows
    staging_batch = []   # pending patric_role_staging rows

    rows_read        = 0
    roles_inserted   = 0
    subs_inserted    = 0
    rs_inserted      = 0
    staging_inserted = 0

    def _flush(start_new_tx=True):
        """Flush pending batches, commit, and optionally start a new transaction."""
        nonlocal rs_inserted, staging_inserted
        if rs_batch:
            conn.executemany(
                "INSERT OR IGNORE INTO role_subsystems(role_id, subsystem_id)"
                " VALUES (?,?)",
                rs_batch,
            )
            rs_inserted += len(rs_batch)
            rs_batch.clear()
        if staging_batch:
            conn.executemany(
                "INSERT INTO patric_role_staging(patric_id, role_id) VALUES (?,?)",
                staging_batch,
            )
            staging_inserted += len(staging_batch)
            staging_batch.clear()
        conn.execute("COMMIT")
        if start_new_tx:
            conn.execute("BEGIN")

    conn.execute("BEGIN")

    with _open_file(subsystem_path) as fh:
        reader = csv.reader(fh, delimiter='\t')
        for row in reader:
            if len(row) < 2:
                continue   # skip malformed rows

            patric_id    = row[0].strip()
            role_name    = row[1].strip()
            superclass   = row[2].strip() if len(row) > 2 else ''
            cls          = row[3].strip() if len(row) > 3 else ''
            subclass     = row[4].strip() if len(row) > 4 else ''
            subsystem_nm = row[5].strip() if len(row) > 5 else ''

            # Required fields
            if not patric_id or not role_name:
                continue

            rows_read += 1

            # ---- role -------------------------------------------------------
            # INSERT OR IGNORE ensures uniqueness; cache avoids re-querying.
            if role_name not in role_cache:
                conn.execute(
                    "INSERT OR IGNORE INTO roles(role_name) VALUES (?)",
                    (role_name,),
                )
                roles_inserted += 1
                role_id = conn.execute(
                    "SELECT role_id FROM roles WHERE role_name = ?",
                    (role_name,),
                ).fetchone()[0]
                role_cache[role_name] = role_id

            role_id = role_cache[role_name]

            # ---- subsystem --------------------------------------------------
            if subsystem_nm:
                if subsystem_nm not in subsystem_cache:
                    conn.execute(
                        "INSERT OR IGNORE INTO subsystems"
                        "(superclass, class, subclass, subsystem_name)"
                        " VALUES (?,?,?,?)",
                        (superclass or None, cls or None,
                         subclass or None, subsystem_nm),
                    )
                    subs_inserted += 1
                    subsystem_id = conn.execute(
                        "SELECT subsystem_id FROM subsystems"
                        " WHERE subsystem_name = ?",
                        (subsystem_nm,),
                    ).fetchone()[0]
                    subsystem_cache[subsystem_nm] = subsystem_id
                else:
                    subsystem_id = subsystem_cache[subsystem_nm]

                rs_key = (role_id, subsystem_id)
                if rs_key not in rs_seen:
                    rs_seen.add(rs_key)
                    rs_batch.append(rs_key)

            staging_batch.append((patric_id, role_id))

            # Flush every batch_size rows to keep memory use bounded.
            if rows_read % batch_size == 0:
                _flush()   # COMMIT + BEGIN new transaction
                if verbose and rows_read % PROGRESS_INTERVAL == 0:
                    print(
                        f"  Subsystem rows processed: {rows_read:,}  "
                        f"roles={len(role_cache):,}  "
                        f"subsystems={len(subsystem_cache):,}",
                        file=sys.stderr,
                    )

    # Final flush – commit remaining rows; do NOT start a new transaction so
    # the next function can open its own.
    _flush(start_new_tx=False)

    if verbose:
        print(
            f"Subsystem file loaded: {rows_read:,} rows  "
            f"roles={roles_inserted:,}  subsystems={subs_inserted:,}  "
            f"role_subsystem_links={rs_inserted:,}  staging={staging_inserted:,}",
            file=sys.stderr,
        )

    return {
        'rows_read':       rows_read,
        'roles_inserted':  roles_inserted,
        'subs_inserted':   subs_inserted,
        'rs_inserted':     rs_inserted,
        'staging_inserted': staging_inserted,
    }


# ---------------------------------------------------------------------------
# Step 2 – load linkout file
# ---------------------------------------------------------------------------

def load_linkout(conn, linkout_path, batch_size, verbose):
    """Stream the patric_uniprot_linkout file into ``linkout_staging``.

    Expected columns (TSV, no header):

    ======  ==========
    col 0   uniprot_id
    col 1   patric_id
    ======  ==========
    """
    batch         = []
    rows_read     = 0
    rows_inserted = 0

    conn.execute("BEGIN")

    with _open_file(linkout_path) as fh:
        reader = csv.reader(fh, delimiter='\t')
        for row in reader:
            if len(row) < 2:
                continue

            uniprot_id = row[0].strip()
            patric_id  = row[1].strip()

            if not uniprot_id or not patric_id:
                continue

            # Gracefully skip an optional header row if present.
            if uniprot_id.lower() in ('uniprot_id', 'uniprot') and \
               patric_id.lower()  in ('patric_id',  'patric'):
                continue

            rows_read += 1
            batch.append((uniprot_id, patric_id))

            if rows_read % batch_size == 0:
                conn.executemany(
                    "INSERT INTO linkout_staging(uniprot_id, patric_id) VALUES (?,?)",
                    batch,
                )
                rows_inserted += len(batch)
                batch.clear()
                conn.execute("COMMIT")
                conn.execute("BEGIN")
                if verbose and rows_read % PROGRESS_INTERVAL == 0:
                    print(
                        f"  Linkout rows processed: {rows_read:,}",
                        file=sys.stderr,
                    )

    if batch:
        conn.executemany(
            "INSERT INTO linkout_staging(uniprot_id, patric_id) VALUES (?,?)",
            batch,
        )
        rows_inserted += len(batch)

    conn.execute("COMMIT")

    if verbose:
        print(
            f"Linkout file loaded: {rows_read:,} rows → "
            f"{rows_inserted:,} staging rows",
            file=sys.stderr,
        )

    return {'rows_read': rows_read, 'rows_inserted': rows_inserted}


# ---------------------------------------------------------------------------
# Step 3 – index staging tables, then JOIN into uniprot_patric_role
# ---------------------------------------------------------------------------

def build_uniprot_patric_role(conn, verbose):
    """Create indexes on the staging tables, then populate ``uniprot_patric_role``
    via a single SQL JOIN between the two staging tables.

    Using a SQL-level JOIN is more efficient than iterating in Python because:
    * SQLite executes the join using the index on patric_role_staging.
    * No Python-level dictionary is needed for the patric_id → role_id lookup.
    * The INSERT...SELECT runs in a single transaction with no per-row overhead.
    """
    if verbose:
        print("Creating staging table indexes...", file=sys.stderr)

    # Create staging indexes AFTER staging data is fully loaded so that bulk
    # inserts above were not slowed down by index maintenance.
    for sql in _INDEX_STAGING:
        conn.execute(sql)

    if verbose:
        print(
            "Populating uniprot_patric_role via staging JOIN...",
            file=sys.stderr,
        )

    conn.execute("BEGIN")
    conn.execute(
        "INSERT OR IGNORE INTO uniprot_patric_role(uniprot_id, patric_id, role_id) "
        "SELECT ls.uniprot_id, ls.patric_id, prs.role_id "
        "FROM   linkout_staging      ls "
        "JOIN   patric_role_staging  prs ON ls.patric_id = prs.patric_id"
    )
    conn.execute("COMMIT")

    n = conn.execute("SELECT COUNT(*) FROM uniprot_patric_role").fetchone()[0]

    if verbose:
        print(f"uniprot_patric_role populated: {n:,} rows", file=sys.stderr)

    return n


# ---------------------------------------------------------------------------
# Step 4 – create final indexes and ANALYZE
# ---------------------------------------------------------------------------

def create_final_indexes(conn, verbose):
    """Create indexes on permanent tables AFTER data load, then run ANALYZE.

    All indexes are deferred until this point so that the large bulk inserts
    above did not incur index-maintenance overhead on every row.
    ANALYZE updates the query-planner statistics so that the new indexes are
    used effectively.
    """
    if verbose:
        print("Creating final indexes on permanent tables...", file=sys.stderr)

    for sql in _INDEX_FINAL:
        conn.execute(sql)

    if verbose:
        print("Running ANALYZE...", file=sys.stderr)

    conn.execute("ANALYZE")


# ---------------------------------------------------------------------------
# Step 5 – optionally drop staging tables
# ---------------------------------------------------------------------------

def drop_staging(conn, verbose):
    """Drop staging tables to reclaim disk space."""
    if verbose:
        print("Dropping staging tables...", file=sys.stderr)
    conn.executescript(_DROP_STAGING)


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def print_summary(conn):
    """Print row counts for all permanent tables to stderr."""
    print("Inserted row counts:", file=sys.stderr)
    for table in ('roles', 'subsystems', 'role_subsystems', 'uniprot_patric_role'):
        n = conn.execute(f"SELECT COUNT(*) FROM [{table}]").fetchone()[0]
        print(f"  {table}: {n:,}", file=sys.stderr)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build a normalized SQLite database from the BV-BRC subsystem "
            "intermediate file and the patric_uniprot linkout file."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        '--linkout', '-l',
        required=True,
        help='patric_uniprot_linkout file (plain .tsv or .gz)',
    )
    parser.add_argument(
        '--subsystems', '-s',
        required=True,
        help='subsystem intermediate TSV file (plain or .gz)',
    )
    parser.add_argument(
        '--db', '-d',
        default='subsystems.db',
        help='output SQLite database path',
    )
    parser.add_argument(
        '--batch-size',
        type=int,
        default=DEFAULT_BATCH_SIZE,
        metavar='N',
        help='number of rows per INSERT batch / transaction commit',
    )
    parser.add_argument(
        '--keep-staging',
        action='store_true',
        default=False,
        help=(
            'keep patric_role_staging and linkout_staging tables after '
            'the final tables have been populated (useful for debugging)'
        ),
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='print progress messages to stderr',
    )

    args = parser.parse_args()

    # --- open / create database ---
    if args.verbose:
        print(f"Opening database: {args.db}", file=sys.stderr)
    conn = _setup_db(args.db)

    # --- step 1: subsystem annotations ---
    if args.verbose:
        print(f"Loading subsystem file: {args.subsystems}", file=sys.stderr)
    load_subsystems(conn, args.subsystems, args.batch_size, args.verbose)

    # --- step 2: linkout ---
    if args.verbose:
        print(f"Loading linkout file: {args.linkout}", file=sys.stderr)
    load_linkout(conn, args.linkout, args.batch_size, args.verbose)

    # --- step 3: populate uniprot_patric_role via staging JOIN ---
    build_uniprot_patric_role(conn, args.verbose)

    # --- step 4: indexes + ANALYZE ---
    create_final_indexes(conn, args.verbose)

    # --- step 5: optional staging cleanup ---
    if not args.keep_staging:
        drop_staging(conn, args.verbose)

    # --- summary ---
    print_summary(conn)

    conn.close()

    if args.verbose:
        print("Done.", file=sys.stderr)


if __name__ == '__main__':
    main()


# ---------------------------------------------------------------------------
# Example queries
# ---------------------------------------------------------------------------

# 1. Look up all roles and subsystem annotations for a given UniProt ID:
#
#    SELECT r.role_name,
#           s.superclass, s.class, s.subclass, s.subsystem_name
#    FROM   uniprot_patric_role upr
#    JOIN   roles           r  ON upr.role_id      = r.role_id
#    JOIN   role_subsystems rs ON r.role_id         = rs.role_id
#    JOIN   subsystems      s  ON rs.subsystem_id   = s.subsystem_id
#    WHERE  upr.uniprot_id = 'P12345';


# 2. Find all UniProt IDs associated with a given subsystem_name:
#
#    SELECT DISTINCT upr.uniprot_id
#    FROM   uniprot_patric_role upr
#    JOIN   role_subsystems rs ON upr.role_id      = rs.role_id
#    JOIN   subsystems      s  ON rs.subsystem_id  = s.subsystem_id
#    WHERE  s.subsystem_name = 'Histidine Biosynthesis';


# 3. Find all roles within a superclass / class / subclass:
#
#    SELECT r.role_name, s.subsystem_name
#    FROM   role_subsystems rs
#    JOIN   roles           r ON rs.role_id      = r.role_id
#    JOIN   subsystems      s ON rs.subsystem_id = s.subsystem_id
#    WHERE  s.superclass = 'Amino Acids and Derivatives'
#      AND  s.class      = 'Histidine Metabolism'
#      AND  s.subclass   = 'Histidine Biosynthesis';


# 4. Count UniProt IDs per subsystem (descending):
#
#    SELECT s.subsystem_name,
#           COUNT(DISTINCT upr.uniprot_id) AS uniprot_count
#    FROM   uniprot_patric_role upr
#    JOIN   role_subsystems rs ON upr.role_id      = rs.role_id
#    JOIN   subsystems      s  ON rs.subsystem_id  = s.subsystem_id
#    GROUP  BY s.subsystem_id
#    ORDER  BY uniprot_count DESC;
