#!/usr/bin/env python3
"""Build deterministic chunk and per-sample merge task files."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--split-root", required=True, type=Path)
    parser.add_argument("--chunk-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--chunk-tasks", required=True, type=Path)
    parser.add_argument("--merge-tasks", required=True, type=Path)
    args = parser.parse_args()

    rows: list[tuple[str, Path, Path]] = []
    for sample_dir in sorted(
        path for path in args.split_root.iterdir() if path.is_dir()
    ):
        if not (sample_dir / "DONE").is_file():
            raise ValueError(f"split lacks DONE marker: {sample_dir}")
        validation = json.loads((sample_dir / "validation.json").read_text())
        with (sample_dir / "manifest.tsv").open() as handle:
            manifest = list(csv.DictReader(handle, delimiter="\t"))
        if validation.get("status") != "validated" or len(manifest) != validation.get(
            "chunks"
        ):
            raise ValueError(f"invalid split metadata: {sample_dir}")
        if sum(int(row["records"]) for row in manifest) != validation["input_records"]:
            raise ValueError(f"record total mismatch: {sample_dir}")
        if sum(int(row["bases"]) for row in manifest) != validation["input_bases"]:
            raise ValueError(f"base total mismatch: {sample_dir}")
        for row in manifest:
            fasta = sample_dir / row["path"]
            if not fasta.is_file() or fasta.stat().st_size == 0:
                raise ValueError(f"missing chunk FASTA: {fasta}")
            label = fasta.name.removesuffix(".fasta.gz")
            rows.append(
                (sample_dir.name, fasta, args.chunk_root / sample_dir.name / label)
            )

    args.chunk_tasks.parent.mkdir(parents=True, exist_ok=True)
    with args.chunk_tasks.open("w") as handle:
        for sample, fasta, output in rows:
            handle.write(f"{sample}__{output.name}\t{fasta}\t{output}\n")
    with args.merge_tasks.open("w") as handle:
        for sample in sorted({sample for sample, _, _ in rows}):
            handle.write(
                f"{sample}\t{args.chunk_root / sample}\t"
                f"{args.chunk_root / sample / 'merged'}\t{args.output_root / sample}\n"
            )


if __name__ == "__main__":
    main()
