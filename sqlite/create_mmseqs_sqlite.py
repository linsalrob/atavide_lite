"""
Import data from mmseqs outputs into sqlite database
"""
import gzip
import os
import sys
import argparse
from atavide_lib import colours

__author__ = 'Rob Edwards'

import sqlite3


def import_data(conn, sample_id, file, table, verbose=False):
    """
    Import the top hits file
    """

    if verbose:
        print(f"{colours.GREEN}Inserting {file} into {table} for sample: {sample_id}{colours.ENDC}", file=sys.stderr)

    opener = open
    if file.endswith('.gz'):
        opener = gzip.open

    cursor = conn.cursor()

    with opener(file, 'rt') as f:
        for l in f:
            p = l.rstrip().split("\t")
            p.insert(0, sample_id)
            if table == 'tophit_aln':
                values = "(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
            elif table == 'lca_taxonomy':
                # k__Bacteria;p__;c__;o__;f__;g__;s__
                taxa = p[-1].split(";")
                p = p + taxa
                if len(p) != 17:
                    ps = "|".join(p)
                    ts = "|".join(taxa)
                    print(f"{colours.RED}Only {len(p)} items in {l}\np: |{ps}\nt: |{ts}|{colours.ENDC}", file=sys.stderr)
                    sys.exit(1)
                values = "(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
            elif table == 'tophit_report_subsystems':
                values = "(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ? ,?, ?)"
            else:
                print(f"{colours.PINK}WARNING: Do not have a SQLite table caled {table}{colours.ENDC}", file=sys.stderr)
                return
            insert_str = f"INSERT INTO {table} VALUES {values}"
            cursor.execute(insert_str, p)

    conn.commit()


def choose_mmseqs(conn, sample_id, mmseqs_dir, verbose=False):
    """
    Import the data from this directory. Currently we import:
    *_tophit_aln.gz
    *__lca_taxonomy.tsv.gz
    *_tophit_report_subsystems.gz
    """

    if verbose:
        print(f"{colours.BLUE}Choosing files for {sample_id}{colours.ENDC}", file=sys.stderr)

    for f in os.listdir(mmseqs_dir):
        if f.endswith('tophit_aln.gz'):
            import_data(conn, sample_id, os.path.join(mmseqs_dir, f), 'tophit_aln', verbose)
        if f.endswith('lca_taxonomy.tsv.gz'):
            import_data(conn, sample_id, os.path.join(mmseqs_dir, f), 'lca_taxonomy', verbose)
        if f.endswith('tophit_report_subsystems.gz'):
            import_data(conn, sample_id, os.path.join(mmseqs_dir, f), 'tophit_report_subsystems', verbose)


def get_all_samples_with_metadata(conn, metadata, mmseqsdir, verbose=False):
    """
    Read the mmseqs dir and get all the sample IDs
    """

    for subdir in os.listdir(mmseqsdir):
        for sample in metadata:
            if subdir.startswith(sample):
                choose_mmseqs(conn, metadata[sample], os.path.join(mmseqsdir, subdir), verbose)


def get_all_samples(conn, mmseqsdir, verbose=False):
    """
    Read the mmseqs dir and get all the sample IDs
    """

    for subdir in os.listdir(mmseqsdir):
        choose_mmseqs(conn, subdir, os.path.join(mmseqsdir, subdir), verbose)



def define_tables(conn, verbose=False):
    if verbose:
        print(f"{colours.BLUE}Defining tables{colours.ENDC}", file=sys.stderr)


    """
    LCA Fields:
    NC_001133.9     4932    species Saccharomyces cerevisiae        32      32      30      0.890   131567;2759;33154;4751;451864;4890;716545;147537;4891;4892;4893;4930;4932
    The contig has been assigned at the species level to taxid 4932 (Saccharomyces cerevisiae). 
    The assignment is based on 32 protein fragments, which passed the fast prefilter selection. 
    Of these, all 32 fragments received a label and of these, 30 agree with the label assigned to the contig 
    (i.e., they were assigned to species 4932 or to a strain below it). The fraction of -log(E-value) support of the label is 89%.
    """

    cursor = conn.cursor()

    # Create the tables
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS tophit_aln (
            sample_id TEXT,
            sequence_id TEXT,
            subject TEXT,
            percent_identity FLOAT,
            alignment_length INTEGER,
            gaps INTEGER,
            mismatches INTEGER,
            query_start INTEGER,
            query_end INTEGER,
            subject_start INTEGER,
            subject_end INTEGER,
            evalue REAL,
            bitscore INTEGER
        );
        """
    )
    conn.commit()
    

    """
    1       R100400180029:20221013142947:V350096418:2:156151:1:2/1/1
    2       1783272
    3       clade
    4       Terrabacteria group
    5       1
    6       1
    7       1
    8       1.000
    9       k__Bacteria;p__;c__;o__;f__;g__;s__
    """

    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS lca_taxonomy (
            sample_id TEXT,
            sequence_id TEXT,
            taxid INTEGER,
            taxonomic_rank TEXT,
            scientific_name TEXT,
            num_fragments INTEGER,
            label_fragments INTEGER,
            label_agreements INTEGER,
            label_support FLOAT,
            phylogeny TEXT,
            taxonomy_kingdom TEXT,
            taxonomy_phylum TEXT,
            taxonomy_class TEXT,
            taxonomy_order TEXT,
            taxonomy_family TEXT,
            taxonomy_genus TEXT,
            taxonomy_species TEXT
        );
        """
    )
    conn.commit()
    

    """
    1       UniRef50_A0A0D6QQV9
    2       1
    3       0.340
    4       0.340
    5       0.367
    6       1300915
    7       species
    8       Anaeromyxobacter sp. PSR-1
    9       Heme A synthase, cytochrome oxidase biogenesis protein Cox15-CtaA
    10      Energy
    11      Respiration
    12      Biogenesis of respiratory chain components
    13      Biogenesis of cytochrome c oxidases
    14      Heme A synthase, cytochrome oxidase biogenesis protein Cox15-CtaA
    15      0.3333333333333333
    """


    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS tophit_report_subsystems (
            sample_id TEXT,
            unirefid TEXT,
            number_aligned_sequences INTEGER,
            unique_target_coverage FLOAT,
            total_target_coverage FLOAT,
            sequence_identity FLOAT,
            taxonomy_id INTEGER,
            taxonomic_rank TEXT,
            scientific_name TEXT,
            uniref_function TEXT,
            area_of_metabolims TEXT,
            classification1 TEXT,
            classification2 TEXT,
            subsystem TEXT,
            rast_function TEXT,
            fraction_of_sequences FLOAT              
        );
        """
    )

    conn.commit()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-d', '--database', help='sqlite database name', required=True)
    parser.add_argument('-m', '--mmseqs', help='directory with mmseqs outputs', required=True)
    parser.add_argument('-x', '--metadata', help='metadata file with subdir\tsample')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    try:
        conn = sqlite3.connect(args.database)
    except sqlite3.Error as e:
        print(f"{colours.RED}FATAL: Can't connect to database {args.database}{colours.ENDC}", file=sys.stderr)
        exit(10)

    define_tables(conn, args.verbose)

    if args.metadata:
        metadata = {}
        with open(args.metadata, 'r') as f:
            for l in f:
                p = l.strip().split("\t")
                metadata[p[0]] = p[1]

        get_all_samples_with_metadata(conn, metadata, args.mmseqs, args.verbose)
    else:
        get_all_samples(conn, args.mmseqs, args.verbose)

    conn.close()
