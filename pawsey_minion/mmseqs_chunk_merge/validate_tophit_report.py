#!/usr/bin/env python3
"""Validate a merged MMseqs top-hit report against all chunk reports."""

from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def read_rows(path: Path):
    seen: set[str] = set()
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(
                    f"{path}:{line_number}: expected at least eight columns"
                )
            if fields[0] in seen:
                raise ValueError(f"{path}:{line_number}: duplicate target {fields[0]}")
            seen.add(fields[0])
            yield fields


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--merged-report", required=True, type=Path)
    parser.add_argument("--chunk-report", required=True, type=Path, action="append")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    chunk_counts: dict[str, int] = {}
    chunk_taxonomy: dict[str, tuple[str, str, str]] = {}
    for report in args.chunk_report:
        for row in read_rows(report):
            target, count, taxonomy = row[0], int(row[1]), tuple(row[5:8])
            if count < 0:
                raise ValueError(f"negative chunk count for {target}")
            if target in chunk_taxonomy and chunk_taxonomy[target] != taxonomy:
                raise ValueError(f"inconsistent chunk taxonomy for {target}")
            chunk_taxonomy[target] = taxonomy
            chunk_counts[target] = chunk_counts.get(target, 0) + count

    merged_targets: set[str] = set()
    merged_count_total = 0
    for row in read_rows(args.merged_report):
        target, count = row[0], int(row[1])
        merged_targets.add(target)
        if count < 0:
            raise ValueError(f"negative merged count for {target}")
        if count != chunk_counts.get(target):
            raise ValueError(f"merged count differs for {target}")
        if tuple(row[5:8]) != chunk_taxonomy[target]:
            raise ValueError(f"merged taxonomy differs for {target}")
        metrics = tuple(map(float, row[2:5]))
        if not all(math.isfinite(value) for value in metrics):
            raise ValueError(f"non-finite metric for {target}")
        if not 0 <= metrics[0] <= 1 or metrics[1] < 0 or not 0 <= metrics[2] <= 1:
            raise ValueError(f"invalid metric for {target}")
        merged_count_total += count
    if merged_targets != set(chunk_counts):
        raise ValueError("target set differs between merged and chunk reports")

    args.output.write_text(
        json.dumps(
            {
                "status": "valid",
                "chunk_reports": len(args.chunk_report),
                "merged_targets": len(merged_targets),
                "merged_count_total": merged_count_total,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
