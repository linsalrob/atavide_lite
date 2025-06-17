"""
Generate the output for sankey_matic to make a sankey plot.
"""

import os
import sys
import argparse
import logging
import gzip
import json
import time
from datetime import datetime
import shutil
from atavide_lib import stream_fastq
from concurrent.futures import ThreadPoolExecutor, as_completed

__author__ = 'Rob Edwards'

# serialise the data as we go along
SERIALISE_EVERY = 30 * 60  # 30 minutes
last_serialised = time.time()

def serialise(data, count_data_file):
    """Serialise the data to a file.
    We do this every SERIALISE_EVERY seconds to avoid losing data."""

    # check the file names and reate the backup files
    if count_data_file.endswith(".gz"):
        backup1 = count_data_file.replace(".gz", ".bak1.gz")
        backup2 = count_data_file.replace(".gz", ".bak2.gz")
    else:
        backup1 = count_data_file + ".bak1.gz"
        backup2 = count_data_file + ".bak2.gz"
        count_data_file += ".gz"

    global last_serialised
    now = time.time()

    if now - last_serialised >= SERIALISE_EVERY:
        readable = datetime.fromtimestamp(now).strftime("%Y-%m-%d %H:%M:%S")
        logging.info(f"Serialising data to {count_data_file} at {readable}")
        # Rotate backups
        if os.path.exists(backup1):
            shutil.move(backup1, backup2)
        if os.path.exists(count_data_file):
            shutil.copy(count_data_file, backup1)

        # Write new JSON
        with gzip.open(count_data_file, "wt") as f:
            json.dump(data, f, indent=2)

        last_serialised = now

def deserialise(count_data_file):
    """
    Read the gziped count file in json format and return the current counts
    """
    if not os.path.exists(count_data_file):
        return {}

    if count_data_file.endswith(".gz"):
        with gzip.open(count_data_file, 'rt') as f:
            return json.load(f)
    else:
        with open(count_data_file, 'r') as f:
            return json.load(f)


def count_fastq(fastq_file, logger=None):
    if logger:
        logger.warning(f"Counting sequences in {fastq_file}")
    count = sum(1 for _ in stream_fastq(fastq_file))
    return count

def count_all_for_read(r):
    """
    Given a read, count the number of sequences in the fastq files,
    but we do it in parallel
    """

    return (
        count_fastq(os.path.join(definitions["SOURCE"], r), logger=logging),
        count_fastq(os.path.join("fastq_fastp", r), logger=logging),
        count_fastq(os.path.join(definitions["HOST"], r), logger=logging),
        count_fastq(os.path.join(definitions["HOSTREMOVED"], r), logger=logging)
    )

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
    parser.add_argument('-s', '--six', help='16S counts file (16S_counts.tsv)', default="16S_counts.tsv")
    parser.add_argument('-c', '--countdata', default="count_data.json.gz",
                        help='intermediate count data file and two backups (default: count_data.json.gz)')
    parser.add_argument('-o', '--output', help='output file', default="sankey_plot.txt")
    parser.add_argument('-j', '--json', help='json output', default='sankey_data.json')
    parser.add_argument('--paired', help='Paired end reads', action='store_true')
    parser.add_argument('-l', '--log', help='log file')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    log_level = "INFO"
    if args.verbose:
        log_level = "DEBUG"

    if args.log:
        logging.basicConfig(filename=args.log, level=log_level, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    else:
        logging.basicConfig(stream=sys.stderr, level=log_level, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    reads = set()
    with open(args.reads, 'r') as reader:
        for line in reader:
            if line.startswith("#") or not line.strip():
                continue
            reads.add(line.strip())

    logging.info(f"Read filess: {len(reads)}")

    definitions = read_defintions(args.definitions, logger=logging)

    reverse_reads = set()
    if args.paired and 'FILEEND' in definitions:
        r2_end = definitions["FILEEND"].replace("R1", "R2")
        for r in reads:
            r2_file = r.replace(definitions["FILEEND"], r2_end)
            reverse_reads.add(r2_file)

    # did we have some data already?
    count_data = {}
    if os.path.exists(args.countdata):
        logging.info(f"Deserialising data from {args.countdata}")
        logging.info(f"Before deserialisation, we have {len(reads.union(reverse_reads))} reads to process")
        count_data = deserialise(args.countdata)

    unprocessed_reads = reads.union(reverse_reads) - set(count_data.keys())
    logging.info(f"Unprocessed reads to count: {len(unprocessed_reads)}")
    
    # read the fastq files
    raw_fastq = 0
    trimmed_fastq = 0
    host = 0
    no_host = 0
    with ThreadPoolExecutor() as executor:
        futures = {executor.submit(count_all_for_read, r): r for r in unprocessed_reads}

        for future in as_completed(futures):
            r = futures[future]
            logging.info(f"Read files for {r}")
            c_raw, c_trimmed, c_host, c_no_host = future.result()
            count_data[r] = {'fastq': c_raw, 'trimmed': c_trimmed, 'host': c_host, 'no_host': c_no_host}
            serialise(count_data, args.countdata)
            raw_fastq += c_raw
            trimmed_fastq += c_trimmed
            host += c_host
            no_host += c_no_host

    # read the 16S counts file if it exists
    logging.info("Reading 16S counts file")
    six_counts = {}
    sixteen_s = 0
    if os.path.exists(args.six):
        with open(args.six, 'r') as six_file:
            for line in six_file:
                if line.startswith("#") or not line.strip() or line.startswith("Sample"):
                    continue
                parts = line.strip().split("\t")
                six_counts[parts[0]] = int(parts[1])
                sixteen_s += int(parts[1])

        logging.info(f"16S counts read: {len(six_counts)}")

    # now read the taxonomy outputs
    tax = {'Bacteria': 0, 'Archaea': 0, 'Eukaryota': 0, 'Viruses': 0, 'Multidomain': 0}
    if os.path.exists(os.path.join("taxonomy_summary", "kingdom.raw.tsv.gz")):
        kingdom_file = "kingdom.raw.tsv.gz" # old format
    elif os.path.exists(os.path.join("taxonomy_summary", f"{definitions['SAMPLENAME']}_kingdom.raw.tsv.gz")):
        kingdom_file = f"{definitions['SAMPLENAME']}_kingdom.raw.tsv.gz" # new format
    else:
        print(f"Error: No kingdom taxonomy file found in {os.path.join('taxonomy_summary')}", file=sys.stderr)

    with gzip.open(os.path.join("taxonomy_summary", kingdom_file), 'rt') as tax_file:
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
            elif 'k__' in parts[0]:
                tax['Multidomain'] = sum(int(x) for x in parts[1:])

    lowqual = raw_fastq - trimmed_fastq
    totalsims = sum(tax.values())
    nosims = no_host-totalsims-sixteen_s


    outputs = f"""
fastq [{trimmed_fastq}] fastp
fastq [{lowqual}] low quality
fastp [{host}] {definitions['HOST']}
fastp [{no_host}] not human
not human [{totalsims}] sequence similarity
not human [{sixteen_s}] 16S
not human [{nosims}] unknown
sequence similarity [{tax['Eukaryota']}] Eukaryote
sequence similarity [{tax['Bacteria']}] Bacteria
sequence similarity [{tax['Archaea']}] Archaea
sequence similarity [{tax['Viruses']}] Virus
sequence similarity [{tax['Multidomain']}] Multiclass
"""

    logging.info(f"Outputs:\n{outputs}")

    with open(args.json, 'w') as jout:
        json.dump({
            "project": definitions.get("SAMPLENAME", "unknown"),
            "raw_fastq": raw_fastq,
            "trimmed_fastq": trimmed_fastq,
            "low_quality": lowqual,
            "host": host,
            "no_host": no_host,
            "sixteen_s": sixteen_s,
            "totalsims": totalsims,
            "nosims": nosims,
            "tax": tax,
            "six_counts": six_counts,
            "count_data": count_data
        }, jout, indent=4)
    with open(args.output, 'w') as out:
        print("Please go to https://sankeymatic.com/build/ and paste this data\n\n", file=out)
        print(outputs, file=out)









