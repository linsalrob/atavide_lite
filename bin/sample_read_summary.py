"""Build a per-sample table of raw read counts and metadata."""

import argparse
import csv
import gzip
import logging
import os
import re
import shlex
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import date


PHROG_REPORT_SUFFIX = "_tophit_report_phrog_function.tsv.gz"
METADATA_COLUMNS = (
    "host_species",
    "individual",
    "sampling_date",
    "water_animal_status",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create a sample-level table of raw pipeline read counts"
    )
    parser.add_argument("-r", "--reads", required=True, help="R1 read-name file")
    parser.add_argument(
        "-d", "--definitions", default="DEFINITIONS.sh", help="definitions file"
    )
    parser.add_argument(
        "-m",
        "--metadata",
        help=(
            "optional tab-separated metadata with a sample column and any of: "
            "host_species, individual, sampling_date, water_animal_status"
        ),
    )
    parser.add_argument(
        "-o", "--output", default="sample_read_summary.tsv", help="output TSV"
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=os.cpu_count() or 1,
        help="parallel sample workers [available CPUs]",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser.parse_args()


def read_definitions(path):
    definitions = {}
    with open(path, encoding="utf-8") as input_file:
        for line_number, line in enumerate(input_file, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("export "):
                line = line[7:].strip()
            if "=" not in line:
                logging.warning("Ignoring %s:%d: no '=' found", path, line_number)
                continue
            key, value = line.split("=", 1)
            try:
                parsed = shlex.split(value, comments=True)
            except ValueError as error:
                raise SystemExit(f"FATAL: cannot parse {path}:{line_number}: {error}")
            definitions[key.strip()] = parsed[0] if parsed else ""
    return definitions


def read_r1_names(path):
    names = []
    seen = set()
    with open(path, encoding="utf-8") as input_file:
        for line_number, line in enumerate(input_file, start=1):
            name = line.strip()
            if not name or name.startswith("#"):
                continue
            if name in seen:
                raise SystemExit(f"FATAL: duplicate read name at {path}:{line_number}: {name}")
            seen.add(name)
            names.append(name)
    if not names:
        raise SystemExit(f"FATAL: no read names were found in {path}")
    return names


def sample_name(read_name, file_end):
    if not read_name.endswith(file_end):
        raise ValueError(f"{read_name!r} does not end with FILEEND={file_end!r}")
    return read_name[: -len(file_end)] if file_end else read_name


def reverse_read_name(read_name, file_end):
    if "R1" not in file_end:
        raise ValueError("FILEEND does not contain R1, so the R2 name cannot be inferred")
    r2_end = file_end.replace("R1", "R2", 1)
    return sample_name(read_name, file_end) + r2_end


def count_fastq(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    opener = gzip.open if path.endswith(".gz") else open
    line_count = 0
    with opener(path, "rb") as input_file:
        for block in iter(lambda: input_file.read(1024 * 1024), b""):
            line_count += block.count(b"\n")
    if line_count % 4:
        raise ValueError(f"{path} has {line_count} lines, which is not a multiple of 4")
    return line_count // 4


def count_pair(directory, r1, r2):
    return count_fastq(os.path.join(directory, r1)) + count_fastq(
        os.path.join(directory, r2)
    )


def parse_phrog_report(path):
    total = 0
    phrog = 0
    with gzip.open(path, "rt", encoding="utf-8") as input_file:
        for line_number, line in enumerate(input_file, start=1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"{path}:{line_number} has fewer than 12 columns")
            try:
                count = int(fields[1])
            except ValueError:
                raise ValueError(
                    f"{path}:{line_number} has a non-integer alignment count: "
                    f"{fields[1]!r}"
                ) from None
            if count < 0:
                raise ValueError(f"{path}:{line_number} has a negative alignment count")
            total += count
            if fields[8]:
                phrog += count
    return total, phrog


def taxonomy_group(label):
    label = label.casefold()
    if "eukary" in label:
        return "Eukarya"
    if "bacteria" in label:
        return "Bacteria"
    if "archaea" in label:
        return "Archaea"
    if "virus" in label or "viridae" in label:
        return "Viruses"
    return "Other"


def parse_taxonomy(path):
    counts = {group: 0 for group in ("Eukarya", "Bacteria", "Archaea", "Viruses", "Other")}
    with gzip.open(path, "rt", encoding="utf-8") as input_file:
        for line_number, line in enumerate(input_file, start=1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number} has fewer than 4 columns")
            try:
                count = int(fields[3])
            except ValueError:
                raise ValueError(
                    f"{path}:{line_number} has a non-integer raw count: {fields[3]!r}"
                ) from None
            if count < 0:
                raise ValueError(f"{path}:{line_number} has a negative raw count")
            counts[taxonomy_group(fields[2])] += count
    return counts


def collect_sample(read_name, definitions):
    file_end = definitions["FILEEND"]
    sample = sample_name(read_name, file_end)
    r2 = reverse_read_name(read_name, file_end)
    source = definitions.get("SOURCE", "fastq")
    host_removed = definitions["HOSTREMOVED"]

    total_reads = count_pair(source, read_name, r2)
    passed_qc = count_pair("fastq_fastp", read_name, r2)
    did_not_map_host = count_pair(host_removed, read_name, r2)

    host_directory = definitions.get("HOST")
    if host_directory:
        mapped_host = count_pair(host_directory, read_name, r2)
    elif os.path.normpath(host_removed) == os.path.normpath("fastq_fastp"):
        mapped_host = 0
    else:
        raise ValueError("HOST is not defined and HOSTREMOVED is not fastq_fastp")

    phrog_report = os.path.join("mmseqs", sample, sample + PHROG_REPORT_SUFFIX)
    mmseqs_hits, phrog_hits = parse_phrog_report(phrog_report)
    taxonomy = parse_taxonomy(os.path.join("taxonomy", sample, "kingdom.tsv.gz"))

    return {
        "sample": sample,
        "total_reads": total_reads,
        "passed_qc_reads": passed_qc,
        "host_mapped_reads": mapped_host,
        "host_unmapped_reads": did_not_map_host,
        "mmseqs_top_hit_reads": mmseqs_hits,
        "phrog_matched_reads": phrog_hits,
        "taxonomically_classified_reads": sum(taxonomy.values()),
        "classified_eukarya_reads": taxonomy["Eukarya"],
        "classified_bacteria_reads": taxonomy["Bacteria"],
        "classified_archaea_reads": taxonomy["Archaea"],
        "classified_viruses_reads": taxonomy["Viruses"],
        "classified_other_reads": taxonomy["Other"],
    }


def normalise_header(value):
    return re.sub(r"[^a-z0-9]+", "_", value.strip().casefold()).strip("_")


def read_metadata(path):
    aliases = {
        "sample": ("sample", "sample_id", "sample_name"),
        "host_species": ("host_species", "species"),
        "individual": ("individual", "individual_id", "animal_id", "host_id"),
        "sampling_date": ("sampling_date", "sample_date", "date"),
        "water_animal_status": (
            "water_animal_status",
            "status",
            "sample_type",
            "source_type",
        ),
    }
    metadata = {}
    with open(path, newline="", encoding="utf-8") as input_file:
        reader = csv.DictReader(input_file, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit(f"FATAL: {path} has no header")
        original_headers = {normalise_header(name): name for name in reader.fieldnames}
        columns = {}
        for destination, candidates in aliases.items():
            columns[destination] = next(
                (original_headers[name] for name in candidates if name in original_headers),
                None,
            )
        if columns["sample"] is None:
            raise SystemExit(f"FATAL: {path} has no sample column")

        for line_number, row in enumerate(reader, start=2):
            sample = (row.get(columns["sample"]) or "").strip()
            if not sample:
                logging.warning("Ignoring %s:%d: empty sample", path, line_number)
                continue
            if sample in metadata:
                raise SystemExit(f"FATAL: duplicate metadata for {sample} in {path}")
            metadata[sample] = {
                field: (row.get(columns[field]) or "").strip() if columns[field] else ""
                for field in METADATA_COLUMNS
            }
    return metadata


def infer_metadata(sample, definitions):
    words = [word for word in re.split(r"[^A-Za-z0-9]+", sample.casefold()) if word]
    joined = "_".join(words)
    inferred = {field: "" for field in METADATA_COLUMNS}

    species_patterns = (
        (r"(?:^|_)homo_sapiens(?:_|$)|(?:^|_)human(?:_|$)", "Homo sapiens"),
        (r"(?:^|_)mus_musculus(?:_|$)|(?:^|_)mouse(?:_|$)", "Mus musculus"),
        (r"(?:^|_)shark(?:_|$)", "shark"),
    )
    for pattern, species in species_patterns:
        if re.search(pattern, joined):
            inferred["host_species"] = species
            break
    if not inferred["host_species"]:
        host = definitions.get("HOST", "")
        if host.casefold() in {"human", "mouse", "shark"}:
            inferred["host_species"] = {
                "human": "Homo sapiens",
                "mouse": "Mus musculus",
                "shark": "shark",
            }[host.casefold()]

    statuses = set()
    if any(word in {"water", "seawater", "freshwater"} for word in words):
        statuses.add("water")
    if any(
        word in {"animal", "host", "tissue", "faecal", "fecal", "stool"}
        for word in words
    ):
        statuses.add("animal")
    if len(statuses) == 1:
        inferred["water_animal_status"] = statuses.pop()

    individual_match = re.search(
        r"(?:^|_)(?:individual|ind|animal)[_-]?([a-z0-9]+)(?:_|$)", joined
    )
    if individual_match:
        inferred["individual"] = individual_match.group(1)

    date_match = re.search(
        r"(?<!\d)(20\d{2})[-_]?([01]\d)[-_]?([0-3]\d)(?!\d)", sample
    )
    if date_match:
        try:
            inferred["sampling_date"] = date(
                *(int(part) for part in date_match.groups())
            ).isoformat()
        except ValueError:
            pass
    return inferred


def main():
    args = parse_args()
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s: %(message)s",
    )
    if args.threads < 1:
        raise SystemExit("FATAL: --threads must be at least 1")

    definitions = read_definitions(args.definitions)
    for required in ("FILEEND", "HOSTREMOVED"):
        if not definitions.get(required):
            raise SystemExit(f"FATAL: {required} is not defined in {args.definitions}")
    read_names = read_r1_names(args.reads)
    supplied_metadata = read_metadata(args.metadata) if args.metadata else {}

    results = {}
    workers = min(args.threads, len(read_names))
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(collect_sample, read_name, definitions): read_name
            for read_name in read_names
        }
        for future in as_completed(futures):
            read_name = futures[future]
            try:
                result = future.result()
            except Exception as error:
                raise SystemExit(f"FATAL: failed to process {read_name}: {error}") from error
            results[result["sample"]] = result
            logging.info("Counted %s", result["sample"])

    fieldnames = [
        "sample",
        "total_reads",
        "passed_qc_reads",
        "host_mapped_reads",
        "host_unmapped_reads",
        "mmseqs_top_hit_reads",
        "phrog_matched_reads",
        "taxonomically_classified_reads",
        "classified_eukarya_reads",
        "classified_bacteria_reads",
        "classified_archaea_reads",
        "classified_viruses_reads",
        "classified_other_reads",
        *METADATA_COLUMNS,
    ]
    with open(args.output, "w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for sample in sorted(results):
            row = results[sample]
            metadata = infer_metadata(sample, definitions)
            for field, value in supplied_metadata.get(sample, {}).items():
                if value:
                    metadata[field] = value
            row.update(metadata)
            writer.writerow(row)


if __name__ == "__main__":
    main()
