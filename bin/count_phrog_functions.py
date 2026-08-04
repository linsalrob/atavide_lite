"""Summarise PHROG-annotated MMseqs top-hit reports across samples."""

import argparse
import csv
import gzip
import math
import os
import re
import sys
from collections import defaultdict


REPORT_SUFFIX = "_tophit_report_phrog_function.tsv.gz"
MARKER_ORDER = (
    "all_phrog_hits",
    "known_function_phrog_hits",
    "unknown_function_phrog_hits",
    "integrase",
    "excisionase",
    "integrase_plus_excisionase",
    "lysogeny_associated_regulators",
    "major_capsid",
    "portal",
    "terminase",
    "tail",
    "holin",
    "endolysin",
    "holin_plus_endolysin",
)

# Marker rules are intentionally conservative and centralized here for review.
EXCISIONASE_EXACT = {
    "excisionase",
    "cox-like excisionase and repressor",
    "excisionase and transcriptional regulator",
    "recombination directionality factor",
    "site-specific recombination directionality factor rdf",
}
LYSOGENY_REGULATOR_EXACT = {
    "ci-like repressor",
    "cro-like repressor",
    "cox-like excisionase and repressor",
    "arc-like repressor",
    "arc-like transcriptional regulator",
    "copg-like transcriptional repressor",
    "transcriptional repressor",
    "anti-repressor",
    "anti-repressor ant",
    "cii-like regulator",
    "cii-like transcriptional activator",
    "ciii anti-termination",
    "immunity protein",
    "immunity to superinfection",
    "excisionase and transcriptional regulator",
}
MAJOR_CAPSID_EXACT = {
    "major capsid protein",
    "major head protein",
    "major head and protease protein",
    "major coat protein",
}
ENDOLYSIN_EXACT = {
    "endolysin",
    "endolysin; inhibits rna polymerase",
    "endolysin; lysm motif",
    "lysozyme domain-containing protein",
    "amidase",
    "internal virion protein with endolysin domain",
    "baseplate hub and tail lysozyme",
    "baseplate hub subunit and tail lysozyme",
    "tail associated lysin",
    "tail fiber protein/ lysozyme",
    "tail protein with lysin activity",
    "minor tail protein with lysin activity",
}
UNKNOWN_ANNOTATIONS = {
    "",
    "-",
    "na",
    "n/a",
    "none",
    "unknown",
    "unknown function",
    "hypothetical protein",
    "uncharacterized protein",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create raw and counts-per-million PHROG and marker tables"
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
    parser.add_argument(
        "-a",
        "--phrog-annotations",
        required=True,
        help=(
            "PHROG annotation TSV with phrog, annot, and category columns; "
            "plain text or gzip-compressed (for example phrog_annot_v4.tsv.gz)"
        ),
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser.parse_args()


def normalise_text(value):
    return " ".join(value.strip().split()).casefold()


def normalise_phrog_id(value):
    value = value.strip()
    match = re.fullmatch(r"(?:phrog_)?(\d+)", value, flags=re.IGNORECASE)
    if not match:
        raise ValueError(
            f"invalid PHROG identifier {value!r}; expected a number or phrog_<number>"
        )
    return str(int(match.group(1)))


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "r", encoding="utf-8", newline="")


def read_phrog_annotations(path):
    if not os.path.isfile(path):
        raise SystemExit(f"FATAL: PHROG annotation database does not exist: {path}")
    if os.path.getsize(path) == 0:
        raise SystemExit(f"FATAL: PHROG annotation database is empty: {path}")
    try:
        with open_text(path) as input_file:
            reader = csv.DictReader(input_file, delimiter="\t")
            required = {"phrog", "annot", "category"}
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                found = ", ".join(reader.fieldnames or []) or "no header"
                raise SystemExit(
                    f"FATAL: {path} must contain tab-separated columns phrog, annot, "
                    f"and category (found: {found})"
                )
            records = {}
            for line_number, row in enumerate(reader, start=2):
                try:
                    phrog_id = normalise_phrog_id(row["phrog"])
                except (AttributeError, ValueError) as error:
                    raise SystemExit(f"FATAL: {path}:{line_number}: {error}") from None
                annotation = (row["annot"] or "").strip()
                category = (row["category"] or "").strip()
                record = (annotation, category)
                if phrog_id in records and records[phrog_id] != record:
                    raise SystemExit(
                        f"FATAL: {path}:{line_number} has a conflicting record for "
                        f"PHROG {phrog_id}: {records[phrog_id]!r} versus {record!r}"
                    )
                records[phrog_id] = record
    except (OSError, UnicodeError, csv.Error) as error:
        raise SystemExit(f"FATAL: cannot read PHROG annotation database {path}: {error}") from None
    if not records:
        raise SystemExit(f"FATAL: PHROG annotation database contains no records: {path}")
    return records


def classify_phrog(annotation, category):
    """Return {marker: auditable rule} for one authoritative PHROG record."""
    annot = normalise_text(annotation)
    cat = normalise_text(category)
    markers = {}

    if re.search(r"\bintegrase\b", annot):
        markers["integrase"] = "word-bounded annotation match: integrase"
    if annot in EXCISIONASE_EXACT:
        markers["excisionase"] = "curated exact excisionase/RDF annotation"
    if annot in LYSOGENY_REGULATOR_EXACT or re.search(
        r"\b(?:repressor|anti-repressor|antirepressor)\b", annot
    ) or annot in {"immunity protein", "immunity to superinfection"}:
        markers["lysogeny_associated_regulators"] = (
            "curated lysogeny regulator or bounded repressor/immunity annotation"
        )
    if annot in MAJOR_CAPSID_EXACT:
        markers["major_capsid"] = "curated exact major capsid/head/coat annotation"
    if re.search(r"\bportal protein\b", annot):
        markers["portal"] = "word-bounded annotation match: portal protein"
    if re.search(r"\bterminase\b", annot):
        markers["terminase"] = "word-bounded annotation match: terminase"
    if cat == "tail":
        markers["tail"] = "exact PHROG category match: tail"
    if re.search(r"\bholin\b", annot):
        markers["holin"] = "word-bounded annotation match: holin"
    tail_or_virion_lysozyme = re.search(
        r"\b(?:virion|tail|baseplate)\b.*\blysozyme\b"
        r"|\blysozyme\b.*\b(?:virion|tail|baseplate)\b",
        annot,
    )
    if annot in ENDOLYSIN_EXACT or re.search(r"\bendolysin\b", annot) or tail_or_virion_lysozyme:
        markers["endolysin"] = (
            "curated muralytic annotation, bounded endolysin, or tail/virion "
            "lysozyme match"
        )

    if "integrase" in markers or "excisionase" in markers:
        markers["integrase_plus_excisionase"] = "union of integrase and excisionase"
    if "holin" in markers or "endolysin" in markers:
        markers["holin_plus_endolysin"] = "union of holin and endolysin"
    return markers


def is_known_function(annotation, category):
    return (
        normalise_text(category) != "unknown function"
        and normalise_text(annotation) not in UNKNOWN_ANNOTATIONS
    )


def add_count(counts, sample, feature, value):
    if feature:
        counts[feature][sample] += value


def format_number(value):
    return f"{value:.12g}"


def write_matrix(path, features, samples, denominator=None, ordered_features=None):
    with open(path, "w", encoding="utf-8") as output:
        output.write("\t" + "\t".join(samples) + "\n")
        feature_names = ordered_features if ordered_features is not None else sorted(features)
        for feature in feature_names:
            output.write(feature)
            for sample in samples:
                value = features[feature].get(sample, 0.0)
                if denominator is not None:
                    total = denominator.get(sample, 0.0)
                    value = value * 1_000_000 / total if total > 0 else 0.0
                output.write("\t" + format_number(value))
            output.write("\n")


def validate_report_metadata(report, line_number, report_value, database_value, field):
    if report_value.strip() and normalise_text(report_value) != normalise_text(database_value):
        raise SystemExit(
            f"FATAL: {report}:{line_number} {field} {report_value!r} conflicts "
            f"with annotation database value {database_value!r}"
        )


def process_report_line(
    line,
    line_number,
    report,
    sample,
    annotation_path,
    annotation_database,
    total_hits,
    phrog_hits,
    phrog_counts,
    annotation_counts,
    category_counts,
    marker_counts,
    metadata,
    observed_database_records,
):
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) < 2:
        raise SystemExit(
            f"FATAL: {report}:{line_number} has {len(fields)} columns; expected at least 2"
        )
    try:
        count = float(fields[1])
    except ValueError:
        raise SystemExit(
            f"FATAL: {report}:{line_number} has a non-numeric alignment count: "
            f"{fields[1]!r}"
        ) from None
    if not math.isfinite(count) or count < 0:
        raise SystemExit(
            f"FATAL: {report}:{line_number} has an invalid alignment count: "
            f"{fields[1]!r}"
        )

    total_hits[sample] += count
    if len(fields) < 12:
        return
    phrog_id, color, report_annotation, report_category = fields[8:12]
    if not phrog_id:
        return
    try:
        canonical_id = normalise_phrog_id(phrog_id)
    except ValueError as error:
        raise SystemExit(f"FATAL: {report}:{line_number}: {error}") from None
    if canonical_id not in annotation_database:
        raise SystemExit(
            f"FATAL: {report}:{line_number} PHROG {phrog_id!r} is absent from "
            f"{annotation_path}"
        )
    db_annotation, db_category = annotation_database[canonical_id]
    validate_report_metadata(
        report, line_number, report_annotation, db_annotation, "annotation"
    )
    validate_report_metadata(
        report, line_number, report_category, db_category, "category"
    )

    phrog_hits[sample] += count
    add_count(phrog_counts, sample, phrog_id, count)
    add_count(annotation_counts, sample, report_annotation, count)
    add_count(category_counts, sample, report_category, count)
    marker_counts["all_phrog_hits"][sample] += count
    known_marker = (
        "known_function_phrog_hits"
        if is_known_function(db_annotation, db_category)
        else "unknown_function_phrog_hits"
    )
    marker_counts[known_marker][sample] += count
    for marker in classify_phrog(db_annotation, db_category):
        marker_counts[marker][sample] += count

    current_metadata = (color, report_annotation, report_category)
    if phrog_id in metadata and metadata[phrog_id] != current_metadata:
        raise SystemExit(
            f"FATAL: {report}:{line_number} has conflicting annotations for {phrog_id}"
        )
    metadata[phrog_id] = current_metadata
    observed_database_records[canonical_id] = (db_annotation, db_category)


def main():
    args = parse_args()
    if not os.path.isdir(args.directory):
        raise SystemExit(f"FATAL: {args.directory} is not a directory")
    annotation_database = read_phrog_annotations(args.phrog_annotations)

    phrog_counts = defaultdict(lambda: defaultdict(float))
    annotation_counts = defaultdict(lambda: defaultdict(float))
    category_counts = defaultdict(lambda: defaultdict(float))
    marker_counts = defaultdict(lambda: defaultdict(float))
    metadata = {}
    observed_database_records = {}
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

        try:
            with gzip.open(report, "rt", encoding="utf-8") as input_file:
                for line_number, line in enumerate(input_file, start=1):
                    process_report_line(
                        line,
                        line_number,
                        report,
                        sample,
                        args.phrog_annotations,
                        annotation_database,
                        total_hits,
                        phrog_hits,
                        phrog_counts,
                        annotation_counts,
                        category_counts,
                        marker_counts,
                        metadata,
                        observed_database_records,
                    )
        except (OSError, UnicodeError) as error:
            raise SystemExit(f"FATAL: cannot read MMseqs report {report}: {error}") from None

    if not samples:
        raise SystemExit(f"FATAL: no *{REPORT_SUFFIX} files were found below {args.directory}")

    os.makedirs(args.output, exist_ok=True)
    levels = {"phrog": phrog_counts, "annotation": annotation_counts, "category": category_counts}
    for level, counts in levels.items():
        write_matrix(os.path.join(args.output, f"{args.name}_{level}_raw.tsv"), counts, samples)
        write_matrix(os.path.join(args.output, f"{args.name}_{level}_norm_all.tsv"), counts, samples, total_hits)
        write_matrix(os.path.join(args.output, f"{args.name}_{level}_norm_phrog.tsv"), counts, samples, phrog_hits)

    write_matrix(os.path.join(args.output, f"{args.name}_marker_raw.tsv"), marker_counts, samples, ordered_features=MARKER_ORDER)
    write_matrix(os.path.join(args.output, f"{args.name}_marker_norm_all.tsv"), marker_counts, samples, total_hits, MARKER_ORDER)
    write_matrix(os.path.join(args.output, f"{args.name}_marker_norm_phrog.tsv"), marker_counts, samples, phrog_hits, MARKER_ORDER)

    with open(os.path.join(args.output, f"{args.name}_phrog_metadata.tsv"), "w", encoding="utf-8") as output:
        output.write("phrog_id\tcolor\tannotation\tcategory\n")
        for phrog_id in sorted(metadata):
            output.write(phrog_id + "\t" + "\t".join(metadata[phrog_id]) + "\n")

    with open(os.path.join(args.output, f"{args.name}_phrog_marker_mapping.tsv"), "w", encoding="utf-8", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(("phrog_id", "annotation", "category", "marker_class", "classification_rule"))
        for canonical_id in sorted(observed_database_records, key=int):
            annotation, category = observed_database_records[canonical_id]
            status = "known_function_phrog_hits" if is_known_function(annotation, category) else "unknown_function_phrog_hits"
            relationships = {"all_phrog_hits": "observed MMseqs top hit with a PHROG match", status: "authoritative annotation/category known-function status"}
            relationships.update(classify_phrog(annotation, category))
            for marker in MARKER_ORDER:
                if marker in relationships:
                    writer.writerow((f"phrog_{canonical_id}", annotation, category, marker, relationships[marker]))

    with open(os.path.join(args.output, f"{args.name}_phrog_count_summary.tsv"), "w", encoding="utf-8") as output:
        output.write("sample\ttotal_mmseqs_top_hits\tall_phrog_hits\tknown_function_phrog_hits\tunknown_function_phrog_hits\tphrog_match_percent\n")
        for sample in samples:
            all_hits = marker_counts["all_phrog_hits"].get(sample, 0.0)
            known = marker_counts["known_function_phrog_hits"].get(sample, 0.0)
            unknown = marker_counts["unknown_function_phrog_hits"].get(sample, 0.0)
            if not math.isclose(known + unknown, all_hits, rel_tol=1e-12, abs_tol=1e-9):
                raise SystemExit(f"FATAL: known and unknown PHROG counts do not sum to all PHROG hits for {sample}")
            percent = all_hits * 100 / total_hits[sample] if total_hits[sample] > 0 else 0.0
            output.write(
                sample
                + "\t"
                + "\t".join(
                    format_number(value)
                    for value in (total_hits[sample], all_hits, known, unknown, percent)
                )
                + "\n"
            )

    with open(os.path.join(args.output, "README.md"), "w", encoding="utf-8") as output:
        output.write(f"""# PHROG count tables

These read-level summaries use the authoritative PHROG annotation database
`{args.phrog_annotations}`. The mapping file includes observed PHROGs only, with
one row per PHROG-marker relationship and the rule that assigned it.

The `{args.name}_*_raw.tsv` files contain alignment counts from column 2 of the
MMseqs top-hit reports. The `{args.name}_*_norm_all.tsv` files contain counts per
1,000,000 all MMseqs top-hit reads, including reads without a PHROG match; this is
the primary normalisation. The `{args.name}_*_norm_phrog.tsv` files contain counts
per 1,000,000 PHROG-matched reads.

Existing tables are provided by PHROG ID, report annotation, and report category.

## Marker definitions

* `all_phrog_hits`: every top hit with a PHROG match.
* `known_function_phrog_hits`: PHROG hits whose authoritative category is not
  `unknown function` and whose annotation is not empty or an unknown placeholder.
* `unknown_function_phrog_hits`: the complementary PHROG-matched hits. Known plus
  unknown therefore equals all PHROG hits for each sample.
* `integrase`: a word-bounded `integrase` annotation.
* `excisionase`: curated exact excisionase and recombination-directionality-factor
  annotations, excluding DNA excision-repair helicases.
* `integrase_plus_excisionase`: the set union of the preceding two classes.
* `lysogeny_associated_regulators`: curated phage repression, immunity,
  establishment, and reversal regulators plus bounded repressor/anti-repressor
  annotations; generic transcriptional regulators and activators are excluded.
* `major_capsid`: curated exact major capsid, major head, and major coat annotations;
  minor and spore coat proteins are excluded.
* `portal`: a word-bounded `portal protein` annotation.
* `terminase`: a word-bounded `terminase` annotation, including large and small
  subunits.
* `tail`: the exact PHROG `tail` category. The `connector` category is not merged.
* `holin`: a word-bounded `holin` annotation, including holin/anti-holin.
* `endolysin`: curated endolysin, lysozyme-domain, amidase, and tail/virion muralytic
  annotations; bounded endolysin and tail/virion-associated lysozyme descriptions
  are also included. Hemolysins and spanins without an independent endolysin
  annotation are excluded.
* `holin_plus_endolysin`: the set union of the preceding two classes.

Combined classes count each read once. Individual component classes can overlap;
for example, a Cox-like excisionase/repressor contributes to both excisionase and
lysogeny-associated regulators, and a tail-associated lysin contributes to tail and
endolysin.

Marker definitions are annotation-based read summaries, not a definitive lytic-
versus-lysogenic classification and not a contig-level lifestyle analysis.
Structural and lysis genes occur in both temperate and virulent phages, and absence
of an integrase does not prove that a phage is virulent.

`{args.name}_phrog_metadata.tsv` preserves metadata observed in reports,
`{args.name}_phrog_marker_mapping.tsv` audits marker assignments, and
`{args.name}_phrog_count_summary.tsv` records denominators and PHROG matching QC.
""")


if __name__ == "__main__":
    main()
