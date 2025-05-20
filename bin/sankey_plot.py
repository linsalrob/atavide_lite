"""
Generate the output for sankey_matic to make a sankey plot.
"""

import os
import sys
import argparse
import logging
from atavide_lib import stream_fastq

__author__ = 'Rob Edwards'

def count_fastq(fastq_file, logger=None):
    if logger:
        logger.warn(f"Counting sequences in {fastq_file}")
    count = sum(1 for _ in stream_fastq(fastq_file))
    return count


def read_defintions(definitions_file, logger=None):
    """
    Read the definitions file and return a dictionary of definitions
    """

    definitions = {}
    with open(definitions_file, 'r') as df:
        for line in df:
            if line.startswith("#") or not line.strip():
                continue
            line = line.replace("export ", "")
            parts = line.strip().split("=")
            if len(parts) != 2:
                print(f"Error: {line.strip()} in {definitions_file} does not have xx=yy parts", file=sys.stderr)
                continue
            definitions[parts[0].strip()] = parts[1].strip()

    logger.debug(f"Read definitions: {definitions}")
    return definitions

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-r', '--reads', help='reads file (e.g. R1_reads.txt)', required=True)
    parser.add_argument('-d', '--definitions', help='DEFINITIONS.sh file', default="DEFINITIONS.sh")
    parser.add_argument('--paired', help='Paired end reads', action='store_true')
    parser.add_argument('-l', '--log', help='log file', default="log.txt")
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    log_level = "INFO"
    if args.verbose:
        log_level = "DEBUG"

    logging.basicConfig(filename=args.log, level=log_level, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    reads = set()
    with open(args.reads, 'r') as reader:
        for line in reader:
            if line.startswith("#") or not line.strip():
                continue
            reads.add(line.strip())

    logging.info(f"Read filess: {len(reads)}")

    definitions = read_defintions(args.definitions, logger=logging)
    r2_end = None
    if 'FILEEND' in definitions:
        r2_end = definitions["FILEEND"].replace("R1", "R2")

    # read the fastq files
    raw_fastq = 0
    trimmed_fastq = 0
    host = 0
    no_host = 0
    for r in reads:
        logging.info(f"Reading files for {r}")
        raw_fastq += count_fastq(os.path.join(definitions["SOURCE"], r), logger=logging)
        trimmed_fastq += count_fastq(os.path.join("fastq_fastp", r), logger=logging)
        host += count_fastq(os.path.join(definitions["HOST"], r), logger=logging)
        no_host += count_fastq(os.path.join(definitions["NO_HOST"], r), logger=logging)

        if r2_end and args.paired:
            logging.info(f"Reading paired end reads for {r}")
            r2_file = r.replace(definitions["FILEEND"], r2_end)
            raw_fastq += count_fastq(os.path.join(definitions["SOURCE"], r2_file), logger=logging)
            trimmed_fastq += count_fastq(os.path.join("fastq_fastp", r2_file), logger=logging)
            host += count_fastq(os.path.join(definitions["HOST"], r2_file), logger=logging)
            no_host += count_fastq(os.path.join(definitions["HOSTREMOVED"], r2_file), logger=logging)

    # now read the taxonomy outputs
    tax = {'Bacteria': 0, 'Archaea': 0, 'Eukaryota': 0, 'Viruses': 0, 'Multidomain': 0}
    if os.path.exists(os.path.join("taxonomy_summary/kingdom.raw.tsv.gz")):
        with open(os.path.join("taxonomy_summary/kingdom.raw.tsv.gz"), 'r') as tax_file:
            for line in tax_file:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if 'Bacteria' in parts[0]:
                    tax['Bacteria'] = sum(int(x) for x in parts[1:])
                elif 'Archaea' in parts[0]:
                    tax['Archaea'] = sum(int(x) for x in parts[1:])
                elif 'Eukaryota' in parts[0]:
                    tax['Eukaryota'] = sum(int(x) for x in parts[1:])
                elif 'Virus' in parts[0]:
                    tax['Viruses'] = sum(int(x) for x in parts[1:])
                else:
                    tax['Multidomain'] = sum(int(x) for x in parts[1:])

    lowqual = raw_fastq - trimmed_fastq
    totalsims = sum(tax.values())
    nosims = no_host-totalsims

    print(f"""
fastq [{trimmed_fastq}] fastp
fastq [{lowqual}] low quality
fastp [{host}] {definitions['HOST']}
fastp [{no_host}] not human
not human [{totalsims}] sequence similarity
not human [{nosims}] unknown
sequence similarity [{tax['Eukaryota']}] Eukaryote
sequence similarity [{tax['Bacteria']}] Bacteria
sequence similarity [{tax['Archaea']}] Archaea
sequence similarity [{tax['Viruses']}] Virus
sequence similarity [{tax['Multidomain']}] Multiclass
""")
