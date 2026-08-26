#!/usr/bin/env python3
"""Generate reproducible, gated VAMB split proposals without replacing bins.

This is deliberately separate from :mod:`refine_hierarchical_vamb`. A forced
latent-space split is merely a proposal: this program evaluates a small k
neighbourhood across deterministic seeds and writes decision-ready evidence.
No assignments or FASTA bins are changed by this script.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import vamb.parsecontigs
import vamb.vambtools

from refine_hierarchical_vamb import marker_calibrated_k, read_candidates, spherical_kmeans_labels


def ari(first: np.ndarray, second: np.ndarray) -> float:
    """Adjusted Rand index, implemented here to avoid a new sklearn dependency."""
    if len(first) != len(second):
        raise ValueError("ARI inputs differ in length")
    _, left = np.unique(first, return_inverse=True)
    _, right = np.unique(second, return_inverse=True)
    table = np.zeros((left.max() + 1, right.max() + 1), dtype=np.int64)
    np.add.at(table, (left, right), 1)
    choose2 = lambda values: np.sum(values * (values - 1) // 2, dtype=np.float64)
    pairs = choose2(table)
    left_pairs = choose2(table.sum(axis=1))
    right_pairs = choose2(table.sum(axis=0))
    total = len(first) * (len(first) - 1) / 2
    if total == 0:
        return 1.0
    expected = left_pairs * right_pairs / total
    maximum = (left_pairs + right_pairs) / 2
    return 1.0 if maximum == expected else float((pairs - expected) / (maximum - expected))


def root(args: argparse.Namespace) -> Path:
    return args.outdir / args.evaluation_dir


def prepare(args: argparse.Namespace) -> None:
    output = root(args)
    inputs = output / "inputs"
    output.mkdir(exist_ok=False)
    inputs.mkdir()
    candidates = read_candidates(args.qa, args.min_completeness, args.min_contamination)
    selected = {str(row["bin"]) for row in candidates}
    contig_parent: dict[str, str] = {}
    with args.assignments.open() as handle:
        next(handle)
        for line in handle:
            parent, contig = line.rstrip("\n").split("\t", 1)
            if parent in selected:
                contig_parent[contig] = parent
    composition = vamb.parsecontigs.Composition.load(args.vamb_dir / "composition.npz")
    identifiers = composition.metadata.identifiers
    lengths = composition.metadata.lengths
    grouped: dict[str, list[int]] = defaultdict(list)
    for index, identifier in enumerate(identifiers):
        parent = contig_parent.get(str(identifier))
        if parent is not None:
            grouped[parent].append(index)
    latent = vamb.vambtools.read_npz(args.vamb_dir / "latent.npz")
    records: list[dict[str, object]] = []
    tasks: list[dict[str, object]] = []
    for candidate in candidates:
        parent = str(candidate["bin"])
        indices = np.asarray(grouped[parent], dtype=np.int64)
        if len(indices) != sum(1 for value in contig_parent.values() if value == parent):
            raise ValueError(f"{parent}: selected contigs missing from composition metadata")
        k0 = marker_calibrated_k(float(candidate["completeness"]), float(candidate["contamination"]), args.max_proposal_k)
        record: dict[str, object] = {
            "bin": parent, "contigs": int(len(indices)), "bases": int(lengths[indices].sum()),
            "completeness": float(candidate["completeness"]), "contamination": float(candidate["contamination"]),
            "k0": k0, "marker_evidence": "PENDING: CheckM per-marker hit locations not retained by node-local QA",
        }
        if k0 > args.max_evaluated_k:
            record["status"] = "MEGABIN_REVIEW"
            records.append(record)
            continue
        record["status"] = "PROPOSED"
        np.savez_compressed(inputs / f"{parent}.npz", names=identifiers[indices], lengths=lengths[indices], latent=latent[indices])
        for k in range(max(2, k0 - 1), min(args.max_evaluated_k, k0 + 1) + 1):
            for seed_index in range(args.seeds):
                tasks.append({"task_id": len(tasks), "bin": parent, "k": k, "seed_index": seed_index})
        records.append(record)
    with (output / "candidates.json").open("w") as handle:
        json.dump({"qa": str(args.qa), "assignments": str(args.assignments), "records": records}, handle, indent=2)
        handle.write("\n")
    with (output / "tasks.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["task_id", "bin", "k", "seed_index"], delimiter="\t")
        writer.writeheader()
        writer.writerows(tasks)
    print(f"prepared candidates={len(records)} evaluation_tasks={len(tasks)} root={output}")


def split(args: argparse.Namespace) -> None:
    output = root(args)
    with (output / "tasks.tsv").open(newline="") as handle:
        tasks = list(csv.DictReader(handle, delimiter="\t"))
    task = next((row for row in tasks if int(row["task_id"]) == args.task_id), None)
    if task is None:
        raise ValueError(f"Unknown task id {args.task_id}")
    parent, k, seed_index = task["bin"], int(task["k"]), int(task["seed_index"])
    with np.load(output / "inputs" / f"{parent}.npz", allow_pickle=True) as arrays:
        names, lengths, latent = arrays["names"], arrays["lengths"], arrays["latent"]
    labels = spherical_kmeans_labels(latent, k, args.seed + args.task_id, args.iterations, args.batch_size, f"cuda:{args.device}")
    target = output / "splits" / parent / f"k{k}"
    target.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(target / f"seed{seed_index}.npz", names=names, lengths=lengths, labels=labels)
    print(f"{parent} k={k} seed={seed_index} contigs={len(names)}")


def summarize(args: argparse.Namespace) -> None:
    output = root(args)
    with (output / "candidates.json").open() as handle:
        records = json.load(handle)["records"]
    rows: list[dict[str, object]] = []
    for record in records:
        parent = str(record["bin"])
        if record["status"] != "PROPOSED":
            rows.append({
                "bin": parent, "k": "", "median_ari": "", "stability": "NOT_EVALUATED",
                "decision": "MEGABIN_REVIEW", "median_ari_ge5kb": "", "ge5kb_contigs": "",
                "medoid_seed": "", "proposal_dir": "", "reason": "k0 exceeds bounded evaluation range",
            })
            continue
        for k_dir in sorted((output / "splits" / parent).glob("k*")):
            paths = sorted(k_dir.glob("seed*.npz"), key=lambda path: int(path.stem.removeprefix("seed")))
            loaded = [np.load(path) for path in paths]
            arrays = [archive["labels"] for archive in loaded]
            if len(arrays) != args.seeds:
                raise ValueError(f"{parent}/{k_dir.name}: expected {args.seeds} seeds, found {len(arrays)}")
            scores = [ari(first, second) for first, second in itertools.combinations(arrays, 2)]
            median = float(np.median(scores))
            names = loaded[0]["names"]
            lengths = loaded[0]["lengths"]
            long_mask = lengths >= args.long_contig_min_length
            long_count = int(long_mask.sum())
            if long_count >= 2:
                long_scores = [
                    ari(first[long_mask], second[long_mask])
                    for first, second in itertools.combinations(arrays, 2)
                ]
                long_median: float | None = float(np.median(long_scores))
            else:
                long_median = None

            mean_agreement = []
            for index, labels in enumerate(arrays):
                agreements = [ari(labels, other) for other_index, other in enumerate(arrays) if other_index != index]
                mean_agreement.append(float(np.mean(agreements)))
            medoid_index = max(range(len(paths)), key=lambda index: (mean_agreement[index], -index))
            medoid_seed = int(paths[medoid_index].stem.removeprefix("seed"))
            stable = (
                median >= args.min_median_ari
                and long_median is not None
                and long_median >= args.min_median_ari
            )
            proposal_dir = ""
            if stable:
                proposal = output / "proposals" / parent / k_dir.name
                proposal.mkdir(parents=True, exist_ok=True)
                labels = arrays[medoid_index]
                with (proposal / "candidate_children.tsv").open("w", newline="") as handle:
                    writer = csv.writer(handle, delimiter="\t")
                    writer.writerow(["clustername", "contigname"])
                    for name, label in zip(names, labels, strict=True):
                        writer.writerow([f"{parent}__{k_dir.name}__child{int(label) + 1}", str(name)])
                child_bases = np.bincount(labels, weights=lengths, minlength=int(k_dir.name[1:])).astype(np.int64)
                child_contigs = np.bincount(labels, minlength=int(k_dir.name[1:])).astype(np.int64)
                with (proposal / "proposal.json").open("w") as handle:
                    json.dump({
                        "parent": parent,
                        "k": int(k_dir.name[1:]),
                        "medoid_seed": medoid_seed,
                        "median_ari_all": median,
                        "median_ari_ge5kb": long_median,
                        "ge5kb_contigs": long_count,
                        "child_contigs": child_contigs.tolist(),
                        "child_bases": child_bases.tolist(),
                        "decision": "PENDING_QC",
                    }, handle, indent=2)
                    handle.write("\n")
                proposal_dir = str(proposal)
            reason = (
                "Requires marker-conflict resolution, improved child quality, and complementary QC"
                if stable else
                "Insufficient contigs >=5 kb for the long-contig stability gate"
                if long_median is None else
                "Split is not reproducible across both all-contig and >=5 kb views"
            )
            rows.append({
                "bin": parent,
                "k": int(k_dir.name[1:]),
                "median_ari": f"{median:.6f}",
                "stability": "STABLE" if stable else "UNSTABLE",
                "decision": "PENDING_QC" if stable else "RETAIN_PARENT",
                "median_ari_ge5kb": "" if long_median is None else f"{long_median:.6f}",
                "ge5kb_contigs": long_count,
                "medoid_seed": medoid_seed,
                "proposal_dir": proposal_dir,
                "reason": reason,
            })
            for archive in loaded:
                archive.close()
    with (output / "split_decisions.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "bin", "k", "median_ari", "stability", "decision", "median_ari_ge5kb",
            "ge5kb_contigs", "medoid_seed", "proposal_dir", "reason",
        ], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {output / 'split_decisions.tsv'}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vamb-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--evaluation-dir", default="gated_refinement_evaluation")
    subparsers = parser.add_subparsers(dest="stage", required=True)
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--qa", type=Path, required=True)
    prepare_parser.add_argument("--assignments", type=Path, required=True)
    prepare_parser.add_argument("--min-completeness", type=float, default=90.0)
    prepare_parser.add_argument("--min-contamination", type=float, default=5.0)
    prepare_parser.add_argument("--max-proposal-k", type=int, default=64)
    prepare_parser.add_argument("--max-evaluated-k", type=int, default=8)
    prepare_parser.add_argument("--seeds", type=int, default=5)
    split_parser = subparsers.add_parser("split")
    split_parser.add_argument("--task-id", type=int, required=True)
    split_parser.add_argument("--seed", type=int, default=60572269627332234)
    split_parser.add_argument("--iterations", type=int, default=30)
    split_parser.add_argument("--batch-size", type=int, default=131_072)
    split_parser.add_argument("--device", type=int, default=0)
    summary_parser = subparsers.add_parser("summarize")
    summary_parser.add_argument("--seeds", type=int, default=5)
    summary_parser.add_argument("--min-median-ari", type=float, default=0.9)
    summary_parser.add_argument("--long-contig-min-length", type=int, default=5000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    {"prepare": prepare, "split": split, "summarize": summarize}[args.stage](args)


if __name__ == "__main__":
    main()
