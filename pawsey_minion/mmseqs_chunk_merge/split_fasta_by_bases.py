#!/usr/bin/env python3
"""Split a FASTA(.gz) into validated, record-complete, base-balanced chunks."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path


def open_text(path: Path, mode: str):
    return gzip.open(path, mode) if path.suffix == ".gz" else path.open(mode)


def records(path: Path):
    header = None
    sequence: list[str] = []
    with open_text(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, sequence
                header, sequence = line, []
            elif header is None:
                raise ValueError("sequence encountered before the first FASTA header")
            else:
                sequence.append(line)
    if header is not None:
        yield header, sequence


def read_id(header: str) -> str:
    return header[1:].split(maxsplit=1)[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--chunks", type=int)
    selection.add_argument("--target-bases", type=int)
    args = parser.parse_args()
    chosen = args.chunks if args.chunks is not None else args.target_bases
    if chosen is None or chosen < 1:
        raise ValueError("chunk count or target bases must be at least one")

    total_bases = total_records = 0
    ids: set[str] = set()
    for header, sequence in records(args.input):
        identifier = read_id(header)
        if identifier in ids:
            raise ValueError(f"duplicate input FASTA identifier: {identifier}")
        ids.add(identifier)
        total_records += 1
        total_bases += sum(len(line.strip()) for line in sequence)
    chunks = args.chunks or max(1, -(-total_bases // args.target_bases))
    if total_records < chunks:
        raise ValueError("cannot create more chunks than FASTA records")

    args.output_dir.mkdir(parents=True, exist_ok=False)
    target = total_bases / chunks
    handles = []
    stats = [
        {"chunk": i + 1, "records": 0, "bases": 0, "ids_sha256": hashlib.sha256()}
        for i in range(chunks)
    ]
    try:
        for i in range(chunks):
            handles.append(
                gzip.open(args.output_dir / f"chunk_{i + 1:03d}.fasta.gz", "wt")
            )
        current = assigned_records = assigned_bases = 0
        for header, sequence in records(args.input):
            bases = sum(len(line.strip()) for line in sequence)
            if (
                current < chunks - 1
                and stats[current]["records"] > 0
                and stats[current]["bases"] + bases > target
                and total_records - assigned_records >= chunks - current - 1
            ):
                current += 1
            handles[current].write(header)
            handles[current].writelines(sequence)
            stats[current]["records"] += 1
            stats[current]["bases"] += bases
            stats[current]["ids_sha256"].update((read_id(header) + "\n").encode())
            assigned_records += 1
            assigned_bases += bases
    finally:
        for handle in handles:
            handle.close()

    emitted_ids: set[str] = set()
    emitted_records = emitted_bases = 0
    for item in stats:
        chunk_path = args.output_dir / f"chunk_{item['chunk']:03d}.fasta.gz"
        for header, sequence in records(chunk_path):
            identifier = read_id(header)
            if identifier in emitted_ids:
                raise AssertionError(
                    f"identifier appears in multiple chunks: {identifier}"
                )
            emitted_ids.add(identifier)
            emitted_records += 1
            emitted_bases += sum(len(line.strip()) for line in sequence)
    if (
        emitted_ids != ids
        or emitted_records != total_records
        or emitted_bases != total_bases
    ):
        raise AssertionError("emitted chunks differ from source FASTA")

    with (args.output_dir / "manifest.tsv").open("w") as handle:
        handle.write("chunk\tpath\trecords\tbases\tids_sha256\n")
        for item in stats:
            handle.write(
                f"{item['chunk']}\tchunk_{item['chunk']:03d}.fasta.gz\t"
                f"{item['records']}\t{item['bases']}\t{item['ids_sha256'].hexdigest()}\n"
            )
    (args.output_dir / "validation.json").write_text(
        json.dumps(
            {
                "input": str(args.input),
                "input_records": total_records,
                "input_bases": total_bases,
                "unique_read_ids": len(ids),
                "chunks": chunks,
                "chunk_records": [item["records"] for item in stats],
                "chunk_bases": [item["bases"] for item in stats],
                "status": "validated",
            },
            indent=2,
        )
        + "\n"
    )
    (args.output_dir / "DONE").write_text(
        "validated record-complete base-balanced split\n"
    )


if __name__ == "__main__":
    main()
