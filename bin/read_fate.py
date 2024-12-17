"""
Figure out the read fate for a set of reads. This is a simple script that will look at the files in a directory and identify
where they went.
"""

import os
import sys
import argparse
import gzip
from atavide_lib import colors, stream_fastq
__author__ = 'Rob Edwards'


def fq_ids(fastq_file):
    """
    Read fqfile and return a set of ids
    """

    ids = set()
    for seqid, header, seq, qual in stream_fastq(fastq_file):
        ids.add(seqid)
    return ids


def sequence_names(read_id_file='R1_reads.txt'):
    """
    Read the read id file and return a set of read filenames
    """

    ids = set()
    with open(read_id_file, 'r') as f:
        for l in f:
            ids.add(l.strip())
    return ids

def read_definitions(defintions_file = "DEFINITIIONS.sh"):
    """
    Read the definitions file and return a dictionary of the definitions
    """

    definitions = {}
    with open(defintions_file, 'r') as f:
        for l in f:
            if l.startswith("export"):
                p = l.replace("export ", "").strip().split("=")
                definitions[p[0]] = p[1].replace('"', '')
    return definitions

def read_tsv(tsvfile, column=0):
    """
    Read a tsv file and return a set of values from the column
    """

    opener = open
    if tsvfile.endswith(".gz"):
        opener = gzip.open
    values = set()
    with opener(tsvfile, 'rt') as f:
        for l in f:
            p = l.strip().split("\t")
            values.add(p[column].replace('/1', '').replace('/2', ''))
    return values


def write_unknown_reads(fastq_file, unknown_reads, outputfile, verbose=False):
    """
    Write the unknown reads to the output file
    """

    if verbose:
        print(f"{colors.GREEN}Reading {fastq_file} and writing unknown reads to {outputfile}{colors.ENDC}")
        print(f"\t{colors.PINK}First two ids :", list(unknown_reads)[0:2], "{colors.ENDC}", file=sys.stderr)
    with open(outputfile, 'a') as out:
        for seqid, header, seq, qual in stream_fastq(fastq_file):
            if seqid in unknown_reads:
                out.write(f"@{header}\n{seq}\n+\n{qual}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-r', '--reads', help='reads file', required=True)
    parser.add_argument('-d', '--definitions', help='definitions file', default='DEFINITIONS.sh')
    parser.add_argument('-u', '--unknown', help='output unknown reads to this file')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()
    
    # step one, read our reads file and get the list of sample names
    full_names = sequence_names(args.reads)

    if args.verbose:
        print(f"{colors.GREEN}Read {len(full_names)} read names{colors.ENDC}")

    # read the definitions file
    defs = read_definitions(args.definitions)
    if args.verbose:
        print(f"{colors.GREEN}Read {len(defs)} definitions{colors.ENDC}")

    short_names = {x: x.replace(defs['FILEEND'], "") for x in full_names}

    # read the initial fastq files
    if args.verbose:
        print(f"{colors.GREEN}Reading fastq files{colors.ENDC}", file=sys.stderr)
    fqfiles = {}
    for n in full_names:
        fqfiles[n] = fq_ids(os.path.join("fastq", n))

    # read the post-qc fastq files
    if args.verbose:
        print(f"{colors.GREEN}Reading post qc fastq files{colors.ENDC}", file=sys.stderr)
    fastp = {}
    for n in full_names:
        fastp[n] = fq_ids(os.path.join("fastq_fastp", n))
    
    # read the human files
    if args.verbose:
        print(f"{colors.GREEN}Reading human files{colors.ENDC}", file=sys.stderr)
    human = {}
    for n in full_names:
        human[n] = fq_ids(os.path.join("human", n))

    # reads the not human files
    if args.verbose:
        print(f"{colors.GREEN}Reading not human files{colors.ENDC}", file=sys.stderr)
    not_human = {}
    for n in full_names:
        not_human[n] = fq_ids(os.path.join("no_human", n))

    # read the mmseqs files
    if args.verbose:
        print(f"{colors.GREEN}Reading mmseqs files{colors.ENDC}", file=sys.stderr)
    mmseqs = {}
    for kn in short_names:
        aln = os.path.join('mmseqs', short_names[kn], f"{short_names[kn]}_tophit_aln.gz")
        if not os.path.exists(aln):
            print(f"{colors.RED}Missing {aln}{colors.ENDC}")
            continue

        mmseqs[kn] = read_tsv(aln, column=0)
        print(f"{colors.GREEN}Read {len(mmseqs[kn])} mmseqs hits for {short_names[kn]}{colors.ENDC}")

    # now we need to figure out where the reads are
    with open("read_fate.tsv", 'w') as out:
        print(f"Name\tTotal\tfastp\tHuman\tNot Human\tMMseqs hits\tUnknown", file=out)
        for n in full_names:
            unknown = len(not_human[n]) - len(mmseqs[n])
            print(f"{short_names[n]}\t{len(fqfiles[n])}\t{len(fastp[n])}\t{len(human[n])}\t{len(not_human[n])}\t{len(mmseqs[n])}\t{unknown}", file=out)

    os.makedirs("read_fate", exist_ok=True)
    # print a list of all the reads in each set
    unknown_reads = {}
    for kn in short_names:
        n = short_names[kn]
        unknown_reads[kn] = set()
        with open(os.path.join("read_fate", f"{n}.txt"), 'w') as f:
            for r in fqfiles[kn]:
                if r not in fastp[kn]:
                    f.write(f"{r}\tfailed fastp\n")
                    continue
                if r in human[kn]:
                    f.write(f"{r}\thuman\n")
                elif r in not_human[kn]:
                    if r in mmseqs[kn]:
                        f.write(f"{r}\tuniref\n")
                    else:
                        unknown_reads[kn].add(r)
                        f.write(f"{r}\tunknown\n")
                else:
                    f.write(f"{r}\tnot sure\n")

    if args.unknown:
        if os.path.exists(args.unknown):
            print(f"{colors.RED}Output file {args.unknown} exists. Cowardly not overwriting it, and we would append to it, so please rename it{colors.ENDC}")
            sys.exit(-1)
        if args.verbose:
            print(f"{colors.GREEN}Writing unknown reads to {args.unknown}{colors.ENDC}", file=sys.stderr)
        for kn in unknown_reads:
            fqfile = os.path.join("fastq", kn)
            write_unknown_reads(fqfile, unknown_reads[kn], args.unknown, args.verbose)
