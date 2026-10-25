#!/usr/bin/env python

import sys
import os
import argparse
import gzip
import hashlib

MAPPING_NAME = "sample_mappings.tsv"
GENERIC_STEMS = {"final.contigs", "contigs", "assembly", "scaffolds"}
GENERIC_DIRS = {"output", "assembly", "megahit"}
FASTA_SUFFIXES = (".gz", ".fa", ".fasta", ".fna")

parser = argparse.ArgumentParser(
    description="""Creates the input FASTA file for Vamb.
Input should be one or more FASTA files, each from a sample-specific assembly.
If keepnames is False, resulting FASTA can be binsplit with separator 'C'.

Vamb renames each contig S<n>C<contig>, where <n> is the position of its input file
on the command line (S1 is the first file). Unless --keepnames is given, the S<n>
to sample mapping is written to sample_mappings.tsv next to the output (see --mapping).
Use --mapping-only to recreate that file for an existing output: the inputs must be
given in the same order, and the mapping is checked against the existing contig names.""",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    add_help=False,
)

parser.add_argument("outpath", help="Path to output FASTA file")
parser.add_argument("inpaths", help="Paths to input FASTA file(s)", nargs="+")
parser.add_argument(
    "-m",
    dest="minlength",
    metavar="",
    type=int,
    default=2000,
    help="Discard sequences below this length [2000]",
)
parser.add_argument(
    "--keepnames", action="store_true", help="Do not rename sequences [False]"
)
parser.add_argument("--nozip", action="store_true", help="Do not gzip output [False]")
parser.add_argument(
    "--mapping",
    metavar="",
    help=f"Path for the S<n> to sample mapping [{MAPPING_NAME} next to outpath]",
)
parser.add_argument(
    "--mapping-only",
    action="store_true",
    help="Only write the mapping for an existing outpath, checking it against the contigs [False]",
)


def sample_name(path):
    """Name a sample from its assembly path.

    megahit/<sample>/final.contigs.fa and megahit/<sample>/output/final.contigs.fa both
    give <sample>; a non-generic file name such as <sample>.fa gives <sample>.
    """
    path = os.path.abspath(path)
    stem = os.path.basename(path)
    while stem.endswith(FASTA_SUFFIXES):
        stem = os.path.splitext(stem)[0]
    if stem not in GENERIC_STEMS:
        return stem
    directory = os.path.dirname(path)
    while os.path.basename(directory) in GENERIC_DIRS:
        directory = os.path.dirname(directory)
    return os.path.basename(directory) or stem


def open_text(path):
    with open(path, "rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if compressed else open(path)


def fasta_records(path):
    """Yield (identifier, sequence) pairs, using the first word of each header as Vamb does."""
    identifier, sequence = None, []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if identifier is not None:
                    yield identifier, "".join(sequence)
                identifier, sequence = line[1:].split(maxsplit=1)[0], []
            else:
                sequence.append(line.strip())
    if identifier is not None:
        yield identifier, "".join(sequence)


def first_kept_contig(path, minlength):
    """Return (identifier, sequence md5) of the first contig Vamb would keep, or None."""
    for identifier, sequence in fasta_records(path):
        if len(sequence) >= minlength:
            return identifier, hashlib.md5(sequence.upper().encode()).hexdigest()
    return None


def check_mapping(outpath, inpaths, minlength):
    """Confirm each input's first kept contig is S<n>C<contig> in outpath, with the same sequence.

    Assemblers reuse contig names across samples (every MEGAHIT assembly has a k141_1),
    so a name alone cannot tell a reordered input list apart; the sequence can.
    """
    expected = {}
    for number, path in enumerate(inpaths, start=1):
        contig = first_kept_contig(path, minlength)
        if contig is not None:
            expected[f"S{number}C{contig[0]}"] = (contig[1], path)
    unmatched = dict(expected)
    for identifier, sequence in fasta_records(outpath):
        if identifier in unmatched:
            if hashlib.md5(sequence.upper().encode()).hexdigest() == unmatched[identifier][0]:
                del unmatched[identifier]
            if not unmatched:
                return
    raise ValueError(
        f"These contigs are missing from {outpath} or have a different sequence there, so "
        "the inputs are not in the order (or -m) used to build it:\n"
        + "\n".join(f"  {name} from {path}" for name, (_, path) in unmatched.items())
    )


def write_mapping(mappath, inpaths):
    with open(mappath, "w") as handle:
        handle.write("vamb_sample\tsample\tpath\n")
        for number, path in enumerate(inpaths, start=1):
            handle.write(f"S{number}\t{sample_name(path)}\t{os.path.abspath(path)}\n")


if len(sys.argv) == 1 or sys.argv[1] in ("-h", "--help"):
    parser.print_help()
    sys.exit()

args = parser.parse_args()

# Check inputs
for path in args.inpaths:
    if not os.path.isfile(path):
        raise FileNotFoundError(path)

if args.keepnames and (args.mapping_only or args.mapping):
    parser.error("--keepnames does not rename contigs, so there is no sample mapping")

mappath = args.mapping or os.path.join(os.path.dirname(args.outpath), MAPPING_NAME)

if args.mapping_only:
    if not os.path.isfile(args.outpath):
        raise FileNotFoundError(args.outpath)
    if os.path.exists(mappath):
        raise FileExistsError(mappath)
    check_mapping(args.outpath, args.inpaths, args.minlength)
    write_mapping(mappath, args.inpaths)
    sys.exit()

if os.path.exists(args.outpath):
    raise FileExistsError(args.outpath)

if not args.keepnames and os.path.exists(mappath):
    raise FileExistsError(mappath)

parent = os.path.dirname(args.outpath)
if parent != "" and not os.path.isdir(parent):
    raise NotADirectoryError(
        f'Output file cannot be created: Parent directory "{parent}" is not an existing directory'
    )

import vamb

# Run the code. Compressing DNA is easy, this is not much bigger than level 9, but
# many times faster
filehandle = (
    open(args.outpath, "w")
    if args.nozip
    else gzip.open(args.outpath, "wt", compresslevel=1)
)
vamb.vambtools.concatenate_fasta(
    filehandle, args.inpaths, minlength=args.minlength, rename=(not args.keepnames)
)
filehandle.close()

if not args.keepnames:
    write_mapping(mappath, args.inpaths)
