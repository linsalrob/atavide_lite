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
import math
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


def read_candidates(qa_path: Path, completeness: float, contamination: float) -> list[dict[str, float | str]]:
    with qa_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Bin Id", "Completeness", "Contamination"}
        if not required.issubset(reader.fieldnames or set()):
            raise ValueError(f"{qa_path} is not a CheckM --tab_table QA report")
        candidates = [
            {
                "bin": row["Bin Id"],
                "completeness": float(row["Completeness"]),
                "contamination": float(row["Contamination"]),
            }
            for row in reader
            if float(row["Completeness"]) >= completeness
            and float(row["Contamination"]) >= contamination
        ]
    return sorted(candidates, key=lambda row: str(row["bin"]))


def marker_calibrated_k(completeness: float, contamination: float, max_clusters: int) -> int:
    """Estimate the number of genomes from duplicate-marker burden.

    CheckM contamination approximates excess copies of lineage marker sets.
    Scaling it by observed completeness gives a conservative lower bound on
    distinct genomes in an otherwise near-complete merged bin.
    """
    estimated = 1 + contamination / max(completeness, 1.0)
    return min(max_clusters, max(2, math.ceil(estimated)))


def spherical_kmeans_labels(
    latent: np.ndarray, n_clusters: int, seed: int, iterations: int, batch_size: int, device: str
) -> np.ndarray:
    """Deterministic forced partitioning of one contaminated parent bin."""
    if n_clusters > len(latent):
        raise ValueError("cluster count exceeds input contig count")
    rows = np.asarray(latent, dtype=np.float32).copy()
    norms = np.linalg.norm(rows, axis=1, keepdims=True)
    rows /= np.where(norms == 0, 1.0, norms)
    rng = np.random.default_rng(seed)
    gpu = torch.device(device)
    centres = torch.from_numpy(rows[rng.choice(len(rows), n_clusters, replace=False)]).to(gpu)
    labels = np.empty(len(rows), dtype=np.int32)
    for iteration in range(iterations):
        sums = torch.zeros_like(centres)
        counts = torch.zeros(n_clusters, dtype=torch.int64, device=gpu)
        for start in range(0, len(rows), batch_size):
            end = min(start + batch_size, len(rows))
            batch = torch.from_numpy(rows[start:end]).to(gpu)
            assigned = (batch @ centres.T).argmax(dim=1)
            labels[start:end] = assigned.cpu().numpy()
            sums.index_add_(0, assigned, batch)
            counts += torch.bincount(assigned, minlength=n_clusters)
        nonempty = counts > 0
        centres[nonempty] = sums[nonempty] / counts[nonempty, None]
        centres[nonempty] /= centres[nonempty].norm(dim=1, keepdim=True)
        if (~nonempty).any():
            replacements = rng.choice(len(rows), int((~nonempty).sum().item()), replace=False)
            centres[~nonempty] = torch.from_numpy(rows[replacements]).to(gpu)
        log(f"forced k-means iteration {iteration + 1}/{iterations}; nonempty={int(nonempty.sum())}")
    return labels


def prepare(args: argparse.Namespace) -> None:
    root = args.outdir / args.refinement_dir
    inputs = root / "inputs"
    root.mkdir(exist_ok=True)
    inputs.mkdir(exist_ok=True)
    candidates = read_candidates(args.qa, args.min_completeness, args.min_contamination)
    if not candidates:
        raise ValueError("No bins meet the refinement thresholds")
    selected = {str(candidate["bin"]) for candidate in candidates}
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
    for candidate in candidates:
        parent = str(candidate["bin"])
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
        forced_k = min(len(indices), marker_calibrated_k(
            float(candidate["completeness"]), float(candidate["contamination"]), args.max_clusters
        ))
        records.append({
            "bin": parent, "contigs": int(len(indices)), "bases": int(lengths[indices].sum()),
            "completeness": float(candidate["completeness"]),
            "contamination": float(candidate["contamination"]), "forced_k": int(forced_k),
        })
        log(f"prepared {parent}: {len(indices):,} contigs, {lengths[indices].sum():,} bp, k={forced_k}")
    (root / "candidates.txt").write_text("\n".join(str(row["bin"]) for row in candidates) + "\n")
    (root / "manifest.json").write_text(json.dumps({
        "qa": str(args.qa), "assignments": str(args.assignments),
        "min_completeness": args.min_completeness,
        "min_contamination": args.min_contamination,
        "candidates": records,
    }, indent=2) + "\n")


def split(args: argparse.Namespace) -> None:
    root = args.outdir / args.refinement_dir
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
    manifest = json.loads((root / "manifest.json").read_text())
    record = next(row for row in manifest["candidates"] if row["bin"] == parent)
    labels = spherical_kmeans_labels(
        latent, int(record["forced_k"]), args.seed + args.task_id,
        args.iterations, args.batch_size, f"cuda:{args.device}"
    )
    emitted = 0
    with output.open("w") as handle:
        print("clustername\tcontigname", file=handle)
        for index, label in enumerate(labels):
            print(f"refined_{parent}_k{int(label) + 1}", names[index], sep="\t", file=handle)
            emitted += 1
    if emitted != len(names):
        raise RuntimeError(f"{parent}: emitted {emitted:,} of {len(names):,} input contigs")
    log(f"{parent}: refined {len(names):,} contigs")


def merge(args: argparse.Namespace) -> None:
    root = args.outdir / args.refinement_dir
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
    parser.add_argument("--refinement-dir", default="refinement")
    subparsers = parser.add_subparsers(dest="stage", required=True)
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--qa", type=Path, required=True)
    prepare_parser.add_argument("--assignments", type=Path, required=True)
    prepare_parser.add_argument("--min-completeness", type=float, default=90.0)
    prepare_parser.add_argument("--min-contamination", type=float, default=5.0)
    prepare_parser.add_argument("--max-clusters", type=int, default=64)
    split_parser = subparsers.add_parser("split")
    split_parser.add_argument("--task-id", type=int, required=True)
    split_parser.add_argument("--seed", type=int, default=60572269627332234)
    split_parser.add_argument("--iterations", type=int, default=30)
    split_parser.add_argument("--batch-size", type=int, default=131_072)
    split_parser.add_argument("--device", type=int, default=0)
    split_parser.add_argument("--force", action="store_true")
    merge_parser = subparsers.add_parser("merge")
    merge_parser.add_argument("--assignments", type=Path, required=True)
    merge_parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.stage == "prepare" and (args.min_completeness < 0 or args.min_contamination < 0 or args.max_clusters < 2):
        parser.error("thresholds must be non-negative")
    return args


def main() -> None:
    args = parse_args()
    {"prepare": prepare, "split": split, "merge": merge}[args.stage](args)


if __name__ == "__main__":
    main()
