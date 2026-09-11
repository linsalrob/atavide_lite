#!/usr/bin/env python3
"""Restore target taxonomy after merging MMseqs swapped-alignment databases."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def open_text(path: Path, mode: str):
    return gzip.open(path, mode) if path.suffix == ".gz" else path.open(mode)


def read_rows(path: Path):
    with open_text(path, "rt") as handle:
        for line_no, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"{path}:{line_no}: expected at least eight columns")
            yield fields


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--merged-report", required=True, type=Path)
    parser.add_argument("--chunk-report", required=True, type=Path, action="append")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    taxonomy: dict[str, tuple[str, str, str]] = {}
    for report in args.chunk_report:
        for row in read_rows(report):
            target, fields = row[0], tuple(row[5:8])
            if target in taxonomy and taxonomy[target] != fields:
                raise ValueError(f"inconsistent chunk taxonomy for {target}")
            taxonomy[target] = fields

    seen: set[str] = set()
    with gzip.open(args.output, "wt") as output:
        for row in read_rows(args.merged_report):
            target = row[0]
            if target in seen:
                raise ValueError(f"duplicate merged target: {target}")
            if target not in taxonomy:
                raise ValueError(f"merged target absent from chunk reports: {target}")
            seen.add(target)
            row[5:8] = taxonomy[target]
            output.write("\t".join(row) + "\n")
    if seen != set(taxonomy):
        raise ValueError("target set differs between merged and chunk reports")


if __name__ == "__main__":
    main()
