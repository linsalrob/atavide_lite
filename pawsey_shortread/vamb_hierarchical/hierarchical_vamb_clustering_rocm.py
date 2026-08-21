#!/usr/bin/env python3
"""Hierarchical, GPU-assisted VAMB clustering for very large latent matrices.

This is intentionally separate from ``resume_vamb_clustering_rocm.py``.
It has three resumable stages:

* ``partition`` learns deterministic spherical k-means centroids and writes
  small, self-contained latent/metadata shards;
* ``local`` runs ordinary global-within-shard VAMB clustering and writes both
  memberships and a weighted representative for every local bin;
* ``collect-representatives`` concatenates those representatives for the
  statistically reconciled second-level clustering step.

The resulting local representatives retain the member mapping, so a later
merge never drops contigs.  Sharding is approximate; use the pilot diagnostics
before interpreting final biological bins.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch
import vamb.cluster
import vamb.parsecontigs
import vamb.vambtools


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def normalized_rows(matrix: np.ndarray) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=np.float32).copy()
    zero = (matrix == 0).all(axis=1)
    matrix[zero] = 1.0 / matrix.shape[1]
    matrix /= np.linalg.norm(matrix, axis=1, keepdims=True) * np.sqrt(2.0)
    return matrix


def spherical_kmeans(
    latent: np.ndarray, n_clusters: int, sample_size: int, iterations: int,
    batch_size: int, seed: int, device: str,
) -> np.ndarray:
    """Fit spherical k-means on a deterministic sample using a GPU."""
    rng = np.random.default_rng(seed)
    sample_n = min(sample_size, len(latent))
    sample_indices = rng.choice(len(latent), size=sample_n, replace=False)
    sample = normalized_rows(latent[sample_indices])
    if sample_n < n_clusters:
        raise ValueError("sample size must be at least the requested shard count")

    gpu = torch.device(device)
    centroids = torch.from_numpy(sample[rng.choice(sample_n, n_clusters, replace=False)]).to(gpu)
    for iteration in range(iterations):
        sums = torch.zeros_like(centroids)
        counts = torch.zeros(n_clusters, dtype=torch.int64, device=gpu)
        for start in range(0, sample_n, batch_size):
            rows = torch.from_numpy(sample[start:start + batch_size]).to(gpu)
            labels = (rows @ centroids.T).argmax(dim=1)
            sums.index_add_(0, labels, rows)
            counts += torch.bincount(labels, minlength=n_clusters)
        nonempty = counts > 0
        centroids[nonempty] = sums[nonempty] / counts[nonempty, None]
        centroids[nonempty] /= centroids[nonempty].norm(dim=1, keepdim=True)
        if (~nonempty).any():
            centroids[~nonempty] = torch.from_numpy(
                sample[rng.choice(sample_n, int((~nonempty).sum().item()), replace=False)]
            ).to(gpu)
        log(f"k-means iteration {iteration + 1}/{iterations}; nonempty={int(nonempty.sum())}")
    return centroids.cpu().numpy()


def assign_shards(latent: np.ndarray, centroids: np.ndarray, batch_size: int, device: str) -> np.ndarray:
    gpu = torch.device(device)
    centres = torch.from_numpy(centroids).to(gpu)
    labels = np.empty(len(latent), dtype=np.int32)
    for start in range(0, len(latent), batch_size):
        end = min(start + batch_size, len(latent))
        rows = torch.from_numpy(normalized_rows(latent[start:end])).to(gpu)
        labels[start:end] = (rows @ centres.T).argmax(dim=1).cpu().numpy()
        if start == 0 or end == len(latent) or (start // batch_size) % 20 == 0:
            log(f"assigned {end:,}/{len(latent):,} contigs")
    return labels


def partition(args: argparse.Namespace) -> None:
    vamb_dir = args.vamb_dir
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    latent = vamb.vambtools.read_npz(vamb_dir / "latent.npz")
    composition = vamb.parsecontigs.Composition.load(vamb_dir / "composition.npz")
    ids, lengths = composition.metadata.identifiers, composition.metadata.lengths
    if len(latent) != len(ids):
        raise ValueError("latent/composition row count mismatch")
    log(f"Partitioning {len(latent):,} contigs into {args.shards} latent-space shards")
    centroids = spherical_kmeans(latent, args.shards, args.sample_size, args.iterations,
                                 args.batch_size, args.seed, args.device)
    labels = assign_shards(latent, centroids, args.batch_size, args.device)
    counts = np.bincount(labels, minlength=args.shards)
    manifest = {"shards": args.shards, "seed": args.seed, "counts": counts.tolist()}
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    np.save(outdir / "centroids.npy", centroids)
    for shard in range(args.shards):
        selected = np.flatnonzero(labels == shard)
        shard_dir = outdir / f"shard_{shard:04d}"
        shard_dir.mkdir(exist_ok=True)
        np.savez_compressed(shard_dir / "latent.npz", latent[selected])
        np.savez_compressed(shard_dir / "metadata.npz", identifiers=ids[selected], lengths=lengths[selected])
    log(f"Wrote {args.shards} shards; sizes min/median/max={counts.min():,}/{int(np.median(counts)):,}/{counts.max():,}")


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


def local(args: argparse.Namespace) -> None:
    shard_dir = args.outdir / f"shard_{args.shard:04d}"
    latent_path, metadata_path = shard_dir / "latent.npz", shard_dir / "metadata.npz"
    if not latent_path.exists():
        raise FileNotFoundError(latent_path)
    output = shard_dir / "local_clusters.tsv"
    reps = shard_dir / "representatives.npz"
    if (output.exists() or reps.exists()) and not args.force:
        raise FileExistsError(f"Refusing to overwrite {shard_dir}; use --force")
    latent = vamb.vambtools.read_npz(latent_path)
    with np.load(metadata_path, allow_pickle=True) as metadata:
        identifiers = metadata["identifiers"]
        lengths = metadata["lengths"]
    if len(latent) == 0:
        return
    torch.cuda.set_device(args.device)
    original = vamb.cluster._calc_distances
    vamb.cluster._calc_distances = chunked_distances(args.chunk_size)
    try:
        generator = vamb.cluster.ClusterGenerator(latent, lengths, destroy=True, cuda=True, rng_seed=args.seed)
        names: list[str] = []
        weights: list[int] = []
        vectors: list[np.ndarray] = []
        with output.open("w") as handle:
            print("clustername\tcontigname", file=handle)
            for number, cluster in enumerate(generator, start=1):
                member_indices = cluster.members.astype(np.int64, copy=False)
                cluster_name = f"{args.shard:04d}:{number}"
                for index in member_indices:
                    print(cluster_name, identifiers[index], sep="\t", file=handle)
                # The CPU array was normalized in place by ClusterGenerator.
                centroid = latent[member_indices].mean(axis=0)
                centroid /= np.linalg.norm(centroid) or 1.0
                names.append(cluster_name)
                weights.append(int(lengths[member_indices].sum()))
                vectors.append(centroid.astype(np.float32, copy=False))
        np.savez_compressed(reps, names=np.array(names, dtype=object), weights=np.array(weights, dtype=np.int64), latent=np.stack(vectors))
        log(f"shard {args.shard}: clustered {len(latent):,} contigs into {len(names):,} local bins")
    finally:
        vamb.cluster._calc_distances = original


def collect_representatives(args: argparse.Namespace) -> None:
    all_names, all_weights, all_latent, all_ncontigs = [], [], [], []
    for shard in range(args.shards):
        shard_dir = args.outdir / f"shard_{shard:04d}"
        path = shard_dir / "representatives.npz"
        if not path.exists():
            raise FileNotFoundError(path)
        with np.load(path, allow_pickle=True) as arrays:
            names = arrays["names"]
            all_names.append(names)
            all_weights.append(arrays["weights"])
            all_latent.append(arrays["latent"])
        # Local cluster names are deterministic ``shard:1``, ``shard:2``, ...
        # in precisely the same order as the representatives. Count members
        # streaming from disk, rather than holding the 22-million-row table.
        ncontigs = np.zeros(len(names), dtype=np.int32)
        memberships = shard_dir / "local_clusters.tsv"
        with memberships.open() as handle:
            next(handle)  # clustername / contigname header
            for line in handle:
                cluster_name = line.partition("\t")[0]
                number = int(cluster_name.partition(":")[2])
                ncontigs[number - 1] += 1
        all_ncontigs.append(ncontigs)
    output = args.outdir / "local_representatives.npz"
    if output.exists() and not args.force:
        raise FileExistsError(f"Refusing to overwrite {output}; use --force")
    n_representatives = sum(map(len, all_names))
    np.savez_compressed(
        output,
        names=np.concatenate(all_names),
        weights=np.concatenate(all_weights),
        latent=np.concatenate(all_latent),
        ncontigs=np.concatenate(all_ncontigs),
    )
    log(f"Wrote {output} with {n_representatives:,} local-bin representatives")


def prepare_global(args: argparse.Namespace) -> None:
    source = args.outdir / "local_representatives.npz"
    output_dir = args.outdir / "global_reconciliation"
    output_dir.mkdir(exist_ok=True)
    with np.load(source, allow_pickle=True) as arrays:
        mask = arrays["ncontigs"] >= args.min_contigs
        names = arrays["names"][mask]
        weights = arrays["weights"][mask]
        latent = arrays["latent"][mask]
        ncontigs = arrays["ncontigs"][mask]
    if len(names) == 0:
        raise ValueError("No local representatives meet --min-contigs")
    np.savez_compressed(output_dir / "latent.npz", latent)
    np.savez_compressed(output_dir / "metadata.npz", names=names, weights=weights, ncontigs=ncontigs)
    log(f"Prepared {len(names):,} representatives with >= {args.min_contigs} contigs")


def cluster_global(args: argparse.Namespace) -> None:
    directory = args.outdir / "global_reconciliation"
    output = directory / "global_bin_members.tsv"
    if output.exists() and not args.force:
        raise FileExistsError(f"Refusing to overwrite {output}; use --force")
    latent = vamb.vambtools.read_npz(directory / "latent.npz")
    with np.load(directory / "metadata.npz", allow_pickle=True) as metadata:
        names = metadata["names"]
        weights = metadata["weights"]
    torch.cuda.set_device(args.device)
    original = vamb.cluster._calc_distances
    vamb.cluster._calc_distances = chunked_distances(args.chunk_size)
    try:
        generator = vamb.cluster.ClusterGenerator(latent, weights, destroy=True, cuda=True, rng_seed=args.seed)
        with output.open("w") as handle:
            print("clustername\tlocal_cluster", file=handle)
            for number, cluster in enumerate(generator, start=1):
                final_name = f"hierarchical_{number}"
                for index in cluster.members:
                    print(final_name, names[int(index)], sep="\t", file=handle)
        log(f"Wrote global reconciliation mapping: {output}")
    finally:
        vamb.cluster._calc_distances = original


def expand_global(args: argparse.Namespace) -> None:
    directory = args.outdir / "global_reconciliation"
    mapping_path = directory / "global_bin_members.tsv"
    output = args.outdir / "hierarchical_unsplit.tsv"
    if output.exists() and not args.force:
        raise FileExistsError(f"Refusing to overwrite {output}; use --force")
    mapping: dict[str, str] = {}
    with mapping_path.open() as handle:
        next(handle)
        for line in handle:
            final_name, local_name = line.rstrip("\n").split("\t", 1)
            mapping[local_name] = final_name
    mapped_contigs = 0
    singleton_contigs = 0
    with output.open("w") as out:
        print("clustername\tcontigname", file=out)
        for shard in range(args.shards):
            memberships = args.outdir / f"shard_{shard:04d}" / "local_clusters.tsv"
            with memberships.open() as handle:
                next(handle)
                for line in handle:
                    local_name, contig_name = line.rstrip("\n").split("\t", 1)
                    final_name = mapping.get(local_name)
                    if final_name is None:
                        # Every representative absent from the global mapping
                        # had one contig, so preserving it as an explicit bin
                        # matches VAMB's loner semantics without dropping it.
                        final_name = f"singleton_{local_name}"
                        singleton_contigs += 1
                    else:
                        mapped_contigs += 1
                    print(final_name, contig_name, sep="\t", file=out)
    log(f"Expanded {mapped_contigs:,} reconciled and {singleton_contigs:,} singleton contigs into {output}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vamb-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    subparsers = parser.add_subparsers(dest="stage", required=True)
    partition_parser = subparsers.add_parser("partition")
    partition_parser.add_argument("--shards", type=int, default=256)
    partition_parser.add_argument("--sample-size", type=int, default=500_000)
    partition_parser.add_argument("--iterations", type=int, default=20)
    partition_parser.add_argument("--batch-size", type=int, default=131_072)
    partition_parser.add_argument("--seed", type=int, default=60572269627332234)
    partition_parser.add_argument("--device", default="cuda:0")
    local_parser = subparsers.add_parser("local")
    local_parser.add_argument("--shard", type=int, required=True)
    local_parser.add_argument("--seed", type=int, default=60572269627332234)
    local_parser.add_argument("--chunk-size", type=int, default=250_000)
    local_parser.add_argument("--device", type=int, default=0)
    local_parser.add_argument("--force", action="store_true")
    collect_parser = subparsers.add_parser("collect-representatives")
    collect_parser.add_argument("--shards", type=int, required=True)
    collect_parser.add_argument("--force", action="store_true")
    prepare_parser = subparsers.add_parser("prepare-global")
    prepare_parser.add_argument("--min-contigs", type=int, default=2)
    global_parser = subparsers.add_parser("cluster-global")
    global_parser.add_argument("--seed", type=int, default=60572269627332234)
    global_parser.add_argument("--chunk-size", type=int, default=250_000)
    global_parser.add_argument("--device", type=int, default=0)
    global_parser.add_argument("--force", action="store_true")
    expand_parser = subparsers.add_parser("expand-global")
    expand_parser.add_argument("--shards", type=int, required=True)
    expand_parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.stage == "partition" and args.shards < 2:
        parser.error("--shards must be at least 2")
    return args


if __name__ == "__main__":
    arguments = parse_args()
    {
        "partition": partition,
        "local": local,
        "collect-representatives": collect_representatives,
        "prepare-global": prepare_global,
        "cluster-global": cluster_global,
        "expand-global": expand_global,
    }[arguments.stage](arguments)
