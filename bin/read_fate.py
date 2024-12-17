"""
Figure out the read fate for a set of reads. This is a simple script that will look at the files in a directory and identify
where they went.
"""

import os
import sys
import argparse
import gzip
from roblib import bcolors, stream_fastq, stream_fasta
__author__ = 'Rob Edwards'

def fq_ids(fqfile, verebose=False):
    """
    Read fqfile and return a set of ids
    """

    ids = set()
    for seqid, header, seq, qual in stream_fastq(fqfile):
        ids.add(seqid)
    return ids


def sequence_names(read_id_file='R1_reads.txt', verbose=False):
    """
    Read the read id file and return a set of read filenames
    """

    ids = set()
    with open(read_id_file, 'r') as f:
        for l in f:
            ids.add(l.strip())
    return ids

def read_definitions(defintions_file = "DEFINITIIONS.sh", verbose=False):
    """
    Read the definitions file and return a dictionary of the definitions
    """

    defs = {}
    with open(defintions_file, 'r') as f:
        for l in f:
            if l.startswith("export"):
                p = l.replace("export ", "").strip().split("=")
                defs[p[0]] = p[1].replace('"', '')
    return defs

def read_tsv(tsvfile, column=0, verbose=False):
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


"""
(miniconda3) edwa0468@setonix-05:Kinshasa [1106]$ ls fastq
barcode01.clean.fastq.gz  barcode03.clean.fastq.gz  barcode05.clean.fastq.gz  barcode07.clean.fastq.gz  barcode09.clean.fastq.gz  barcode11.clean.fastq.gz  barcode13.clean.fastq.gz
barcode02.clean.fastq.gz  barcode04.clean.fastq.gz  barcode06.clean.fastq.gz  barcode08.clean.fastq.gz  barcode10.clean.fastq.gz  barcode12.clean.fastq.gz  unclassified.clean.fastq.gz
(miniconda3) edwa0468@setonix-05:Kinshasa [1107]$ ls human/
barcode01.clean.fastq.gz  barcode03.clean.fastq.gz  barcode05.clean.fastq.gz  barcode07.clean.fastq.gz  barcode09.clean.fastq.gz  barcode11.clean.fastq.gz  barcode13.clean.fastq.gz
barcode02.clean.fastq.gz  barcode04.clean.fastq.gz  barcode06.clean.fastq.gz  barcode08.clean.fastq.gz  barcode10.clean.fastq.gz  barcode12.clean.fastq.gz  unclassified.clean.fastq.gz
(miniconda3) edwa0468@setonix-05:Kinshasa [1109]$ ls mmseqs
barcode01  barcode02  barcode03  barcode04  barcode05  barcode06  barcode07  barcode08  barcode09  barcode10  barcode11  barcode12  barcode13  unclassified
(miniconda3) edwa0468@setonix-05:Kinshasa [1110]$ ls mmseqs/barcode01/
barcode01_lca_taxonomy.tsv.gz  barcode01_lca.tsv.gz  barcode01_report.gz  barcode01_taxonomy.tsv.gz  barcode01_tophit_aln.gz  barcode01_tophit_report.gz  barcode01_tophit_report_subsystems_taxa.gz
(miniconda3) edwa0468@setonix-05:Kinshasa [1111]$
"""



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-r', '--reads', help='reads file', required=True)
    parser.add_argument('-d', '--definitions', help='definitions file', default='DEFINITIONS.sh')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()
    
    # step one, read our reads file and get the list of sample names
    full_names = sequence_names(args.reads, args.verbose)

    if args.verbose:
        print(f"{bcolors.GREEN}Read {len(full_names)} read names{bcolors.ENDC}")

    # read the definitions file
    defs = read_definitions(args.definitions, args.verbose)
    if args.verbose:
        print(f"{bcolors.GREEN}Read {len(defs)} definitions{bcolors.ENDC}")

    short_names = {x:x.replace(defs['FILEEND'], "") for x in full_names}

    # read the initial fastq files
    if args.verbose:
        print(f"{bcolors.GREEN}Reading fastq files{bcolors.ENDC}", file=sys.stderr)
    fqfiles = {}
    for n in full_names:
        fqfiles[n] = fq_ids(os.path.join("fastq", n), args.verbose)

    # read the post-qc fastq files
    if args.verbose:
        print(f"{bcolors.GREEN}Reading post qc fastq files{bcolors.ENDC}", file=sys.stderr)
    fastp = {}
    for n in full_names:
        fastp[n] = fq_ids(os.path.join("fastq_fastp", n), args.verbose)
    
    # read the human files
    if args.verbose:
        print(f"{bcolors.GREEN}Reading human files{bcolors.ENDC}", file=sys.stderr)
    human = {}
    for n in full_names:
        human[n] = fq_ids(os.path.join("human", n), args.verbose)

    # reads the not human files
    if args.verbose:
        print(f"{bcolors.GREEN}Reading not human files{bcolors.ENDC}", file=sys.stderr)
    not_human = {}
    for n in full_names:
        not_human[n] = fq_ids(os.path.join("no_human", n), args.verbose)

    # read the mmseqs files
    if args.verbose:
        print(f"{bcolors.GREEN}Reading mmseqs files{bcolors.ENDC}", file=sys.stderr)
    mmseqs = {}
    for kn in short_names:
        aln = os.path.join('mmseqs', short_names[kn], f"{short_names[kn]}_tophit_aln.gz")
        if not os.path.exists(aln):
            print(f"{bcolors.RED}Missing {aln}{bcolors.ENDC}")
            continue

        mmseqs[kn] = read_tsv(aln, column=0, verbose=args.verbose)
        print(f"{bcolors.GREEN}Read {len(mmseqs[kn])} mmseqs hits for {short_names[kn]}{bcolors.ENDC}")


    # now we need to figure out where the reads are

    with open("read_fate.tsv", 'w') as out:
        print(f"Name\tTotal\tfastp\tHuman\tNot Human\tMMseqs hits\tUnknown", file=out)
        for n in full_names:
            unknown = len(not_human[n]) - len(mmseqs[n])
            print(f"{short_names[n]}\t{len(fqfiles[n])}\t{len(fastp[n])}\t{len(human[n])}\t{len(not_human[n])}\t{len(mmseqs[n])}\t{unknown}", file=out)

    os.makedirs("read_fate", exist_ok=True)
    # print a list of all the reads in each set
    for kn in short_names:
        n = short_names[kn]
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
                        f.write(f"{r}\tunknown\n")
                else:
                    f.write(f"{r}\tnot sure\n")


