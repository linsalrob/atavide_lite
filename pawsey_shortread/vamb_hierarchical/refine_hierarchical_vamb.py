#!/usr/bin/env python3
"""CheckM-guided, targeted refinement of contaminated hierarchical VAMB bins.

The workflow deliberately acts only on bins selected from a CheckM tabular QA
report. It reclusters the original VAMB contig embeddings for each selected
bin, then expands the refined memberships back into a complete assignment
table without discarding untreated bins or contigs.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch
import vamb.cluster
import vamb.parsecontigs
import vamb.vambtools


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def chunked_distances(chunk_size: int):
    def calculate(matrix: torch.Tensor, index: int) -> torch.Tensor:
        vector = matrix[index]
        distances = torch.empty(len(matrix), device=matrix.device, dtype=matrix.dtype)
        for start in range(0, len(matrix), chunk_size):
            end = min(start + chunk_size, len(matrix))
            distances[start:end] = 0.5 - torch.mv(matrix[start:end], vector)
        distances[index] = 0.0
        return distances

    return calculate


def read_candidates(qa_path: Path, completeness: float, contamination: float) -> list[str]:
    with qa_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Bin Id", "Completeness", "Contamination"}
        if not required.issubset(reader.fieldnames or set()):
            raise ValueError(f"{qa_path} is not a CheckM --tab_table QA report")
        candidates = [
            row["Bin Id"]
            for row in reader
            if float(row["Completeness"]) >= completeness
            and float(row["Contamination"]) >= contamination
        ]
    return sorted(candidates)


def prepare(args: argparse.Namespace) -> None:
    root = args.outdir / "refinement"
    inputs = root / "inputs"
    root.mkdir(exist_ok=True)
    inputs.mkdir(exist_ok=True)
    candidates = read_candidates(args.qa, args.min_completeness, args.min_contamination)
    if not candidates:
        raise ValueError("No bins meet the refinement thresholds")
    selected = set(candidates)
    contig_parent: dict[str, str] = {}
    with args.assignments.open() as handle:
        next(handle)
        for line in handle:
            bin_name, contig = line.rstrip("\n").split("\t", 1)
            if bin_name in selected:
                if contig in contig_parent:
                    raise ValueError(f"Contig occurs more than once in assignments: {contig}")
                contig_parent[contig] = bin_name
    observed = set(contig_parent.values())
    missing = selected - observed
    if missing:
        raise ValueError(f"QA bins absent from assignments: {', '.join(sorted(missing)[:5])}")

    composition = vamb.parsecontigs.Composition.load(args.vamb_dir / "composition.npz")
    identifiers = composition.metadata.identifiers
    lengths = composition.metadata.lengths
    grouped_indices: dict[str, list[int]] = defaultdict(list)
    for index, identifier in enumerate(identifiers):
        parent = contig_parent.get(str(identifier))
        if parent is not None:
            grouped_indices[parent].append(index)
    if sum(map(len, grouped_indices.values())) != len(contig_parent):
        raise ValueError("Some selected contigs were absent from VAMB composition metadata")

    latent = vamb.vambtools.read_npz(args.vamb_dir / "latent.npz")
    records = []
    for parent in candidates:
        indices = np.asarray(grouped_indices[parent], dtype=np.int64)
        if len(indices) < 2:
            raise ValueError(f"Cannot refine bin with fewer than two contigs: {parent}")
        path = inputs / f"{parent}.npz"
        np.savez_compressed(
            path,
            names=identifiers[indices],
            lengths=lengths[indices],
            latent=latent[indices],
        )
        records.append({"bin": parent, "contigs": int(len(indices)), "bases": int(lengths[indices].sum())})
        log(f"prepared {parent}: {len(indices):,} contigs, {lengths[indices].sum():,} bp")
    (root / "candidates.txt").write_text("\n".join(candidates) + "\n")
    (root / "manifest.json").write_text(json.dumps({
        "qa": str(args.qa), "assignments": str(args.assignments),
        "min_completeness": args.min_completeness,
        "min_contamination": args.min_contamination,
        "candidates": records,
    }, indent=2) + "\n")


def split(args: argparse.Namespace) -> None:
    root = args.outdir / "refinement"
    candidates = (root / "candidates.txt").read_text().splitlines()
    if args.task_id < 0 or args.task_id >= len(candidates):
        raise ValueError(f"task id {args.task_id} outside candidate range 0..{len(candidates) - 1}")
    parent = candidates[args.task_id]
    output_dir = root / "splits"
    output_dir.mkdir(exist_ok=True)
    output = output_dir / f"{parent}.tsv"
    if output.exists() and not args.force:
        raise FileExistsError(f"Refusing to overwrite {output}; use --force")
    with np.load(root / "inputs" / f"{parent}.npz", allow_pickle=True) as arrays:
        names = arrays["names"]
        lengths = arrays["lengths"]
        latent = arrays["latent"]
    torch.cuda.set_device(args.device)
    original = vamb.cluster._calc_distances
    vamb.cluster._calc_distances = chunked_distances(args.chunk_size)
    emitted = 0
    try:
        generator = vamb.cluster.ClusterGenerator(latent, lengths, destroy=True, cuda=True, rng_seed=args.seed)
        with output.open("w") as handle:
            print("clustername\tcontigname", file=handle)
            for number, cluster in enumerate(generator, start=1):
                name = f"refined_{parent}_{number}"
                for index in cluster.members.astype(np.int64, copy=False):
                    print(name, names[index], sep="\t", file=handle)
                    emitted += 1
    finally:
        vamb.cluster._calc_distances = original
    if emitted != len(names):
        raise RuntimeError(f"{parent}: emitted {emitted:,} of {len(names):,} input contigs")
    log(f"{parent}: refined {len(names):,} contigs")


def merge(args: argparse.Namespace) -> None:
    root = args.outdir / "refinement"
    candidates = (root / "candidates.txt").read_text().splitlines()
    replacement: dict[str, str] = {}
    for parent in candidates:
        path = root / "splits" / f"{parent}.tsv"
        if not path.exists():
            raise FileNotFoundError(path)
        with path.open() as handle:
            next(handle)
            for line in handle:
                refined, contig = line.rstrip("\n").split("\t", 1)
                if contig in replacement:
                    raise ValueError(f"Contig appears in more than one refinement: {contig}")
                replacement[contig] = refined
    output = root / "refined_unsplit.tsv"
    if output.exists() and not args.force:
        raise FileExistsError(f"Refusing to overwrite {output}; use --force")
    replaced = 0
    with args.assignments.open() as source, output.open("w") as destination:
        print(next(source).rstrip("\n"), file=destination)
        for line in source:
            parent, contig = line.rstrip("\n").split("\t", 1)
            refined = replacement.get(contig)
            if refined is None:
                print(parent, contig, sep="\t", file=destination)
            else:
                print(refined, contig, sep="\t", file=destination)
                replaced += 1
    if replaced != len(replacement):
        raise RuntimeError(f"Merged {replaced:,} of {len(replacement):,} refined contigs")
    log(f"Wrote {output}; replaced {replaced:,} contigs across {len(candidates):,} parent bins")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vamb-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    subparsers = parser.add_subparsers(dest="stage", required=True)
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--qa", type=Path, required=True)
    prepare_parser.add_argument("--assignments", type=Path, required=True)
    prepare_parser.add_argument("--min-completeness", type=float, default=90.0)
    prepare_parser.add_argument("--min-contamination", type=float, default=5.0)
    split_parser = subparsers.add_parser("split")
    split_parser.add_argument("--task-id", type=int, required=True)
    split_parser.add_argument("--seed", type=int, default=60572269627332234)
    split_parser.add_argument("--chunk-size", type=int, default=250_000)
    split_parser.add_argument("--device", type=int, default=0)
    split_parser.add_argument("--force", action="store_true")
    merge_parser = subparsers.add_parser("merge")
    merge_parser.add_argument("--assignments", type=Path, required=True)
    merge_parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.stage == "prepare" and (args.min_completeness < 0 or args.min_contamination < 0):
        parser.error("thresholds must be non-negative")
    return args


def main() -> None:
    args = parse_args()
    {"prepare": prepare, "split": split, "merge": merge}[args.stage](args)


if __name__ == "__main__":
    main()
