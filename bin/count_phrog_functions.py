"""Summarise PHROG-annotated MMseqs top-hit reports across samples."""

import argparse
import gzip
import math
import os
import sys
from collections import defaultdict


REPORT_SUFFIX = "_tophit_report_phrog_function.tsv.gz"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create raw and counts-per-million PHROG summary tables"
    )
    parser.add_argument(
        "-d",
        "--directory",
        default="mmseqs",
        help="directory containing one MMseqs output directory per sample [mmseqs]",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="phrog_functions",
        help="output directory [phrog_functions]",
    )
    parser.add_argument(
        "-n", "--name", default="atavide", help="prefix for output files [atavide]"
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser.parse_args()


def add_count(counts, sample, feature, value):
    if feature:
        counts[feature][sample] += value


def format_number(value):
    return f"{value:.12g}"


def write_matrix(path, features, samples, denominator=None):
    with open(path, "w", encoding="utf-8") as output:
        output.write("\t" + "\t".join(samples) + "\n")
        for feature in sorted(features):
            output.write(feature)
            for sample in samples:
                value = features[feature].get(sample, 0.0)
                if denominator is not None:
                    total = denominator.get(sample, 0.0)
                    value = value * 1_000_000 / total if total > 0 else 0.0
                output.write("\t" + format_number(value))
            output.write("\n")


def main():
    args = parse_args()
    if not os.path.isdir(args.directory):
        raise SystemExit(f"FATAL: {args.directory} is not a directory")

    phrog_counts = defaultdict(lambda: defaultdict(float))
    annotation_counts = defaultdict(lambda: defaultdict(float))
    category_counts = defaultdict(lambda: defaultdict(float))
    metadata = {}
    total_hits = defaultdict(float)
    phrog_hits = defaultdict(float)
    samples = []

    for sample in sorted(os.listdir(args.directory)):
        sample_directory = os.path.join(args.directory, sample)
        if not os.path.isdir(sample_directory):
            continue

        report = os.path.join(sample_directory, sample + REPORT_SUFFIX)
        if not os.path.isfile(report):
            if args.verbose:
                print(f"Skipping {sample}: {report} does not exist", file=sys.stderr)
            continue

        if args.verbose:
            print(f"Reading {sample} from {report}", file=sys.stderr)
        samples.append(sample)

        with gzip.open(report, "rt", encoding="utf-8") as input_file:
            for line_number, line in enumerate(input_file, start=1):
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 12:
                    raise SystemExit(
                        f"FATAL: {report}:{line_number} has {len(fields)} columns; "
                        "expected 12"
                    )

                try:
                    count = float(fields[1])
                except ValueError:
                    raise SystemExit(
                        f"FATAL: {report}:{line_number} has a non-numeric alignment "
                        f"count: {fields[1]!r}"
                    ) from None
                if not math.isfinite(count) or count < 0:
                    raise SystemExit(
                        f"FATAL: {report}:{line_number} has an invalid alignment "
                        f"count: {fields[1]!r}"
                    )

                total_hits[sample] += count
                phrog_id, color, annotation, category = fields[8:12]
                if not phrog_id:
                    continue

                phrog_hits[sample] += count
                add_count(phrog_counts, sample, phrog_id, count)
                add_count(annotation_counts, sample, annotation, count)
                add_count(category_counts, sample, category, count)

                current_metadata = (color, annotation, category)
                if phrog_id in metadata and metadata[phrog_id] != current_metadata:
                    raise SystemExit(
                        f"FATAL: {report}:{line_number} has conflicting annotations "
                        f"for {phrog_id}"
                    )
                metadata[phrog_id] = current_metadata

    if not samples:
        raise SystemExit(
            f"FATAL: no *{REPORT_SUFFIX} files were found below {args.directory}"
        )

    os.makedirs(args.output, exist_ok=True)
    levels = {
        "phrog": phrog_counts,
        "annotation": annotation_counts,
        "category": category_counts,
    }
    for level, counts in levels.items():
        write_matrix(
            os.path.join(args.output, f"{args.name}_{level}_raw.tsv"),
            counts,
            samples,
        )
        write_matrix(
            os.path.join(args.output, f"{args.name}_{level}_norm_all.tsv"),
            counts,
            samples,
            total_hits,
        )
        write_matrix(
            os.path.join(args.output, f"{args.name}_{level}_norm_phrog.tsv"),
            counts,
            samples,
            phrog_hits,
        )

    with open(
        os.path.join(args.output, f"{args.name}_phrog_metadata.tsv"),
        "w",
        encoding="utf-8",
    ) as output:
        output.write("phrog_id\tcolor\tannotation\tcategory\n")
        for phrog_id in sorted(metadata):
            output.write(phrog_id + "\t" + "\t".join(metadata[phrog_id]) + "\n")

    with open(
        os.path.join(args.output, "README.md"), "w", encoding="utf-8"
    ) as output:
        output.write(
            f"""# PHROG count tables

The `{args.name}_*_raw.tsv` files contain alignment counts from column 2 of the
MMseqs top-hit reports.

The `{args.name}_*_norm_all.tsv` files contain counts per million of all MMseqs
top hits, including hits without a PHROG match.

The `{args.name}_*_norm_phrog.tsv` files contain counts per million of MMseqs
top hits that have a PHROG match.

Tables are provided by PHROG ID, annotation, and category. The
`{args.name}_phrog_metadata.tsv` file maps each observed PHROG ID to its color,
annotation, and category.
"""
        )


if __name__ == "__main__":
    main()
