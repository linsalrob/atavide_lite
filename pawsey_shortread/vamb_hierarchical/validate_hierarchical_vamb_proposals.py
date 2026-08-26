#!/usr/bin/env python3
"""Materialize and assess stable hierarchical VAMB split proposals.

This program never changes canonical assignments. It materializes the medoid
seed selected by ``evaluate_hierarchical_vamb_refinement.py`` and evaluates
parent and child bins with the *same parent-derived CheckM marker set*.
Marker resolution is necessary but not sufficient: no result produced here is
an accepted split without complementary QC.
"""
from __future__ import annotations

import argparse
import ast
import csv
import json
import shutil
from collections import defaultdict
from pathlib import Path

def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def materialize(args: argparse.Namespace) -> None:
    import vamb.vambtools

    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite marker-validation directory: {args.output}")
    bins_dir = args.output / "bins"
    parents_dir = args.output / "parents"
    bins_dir.mkdir(parents=True)
    parents_dir.mkdir()
    pending = [row for row in read_tsv(args.decisions) if row["decision"] == "PENDING_QC"]
    by_parent: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in pending:
        by_parent[row["bin"]].append(row)

    manifest: list[dict[str, object]] = []
    for parent, proposals in sorted(by_parent.items()):
        parent_fasta = args.parent_bins / f"{parent}.fna"
        if not parent_fasta.is_file():
            raise FileNotFoundError(f"Missing canonical parent FASTA: {parent_fasta}")
        parent_ids: set[str] = set()
        with vamb.vambtools.Reader(parent_fasta) as reader:
            for entry in vamb.vambtools.byte_iterfasta(reader, str(parent_fasta)):
                parent_ids.add(entry.identifier)
        original_id = f"{parent}__original"
        original_fasta = bins_dir / f"{original_id}.fna"
        shutil.copyfile(parent_fasta, original_fasta)
        (parents_dir / original_fasta.name).symlink_to(Path("..") / "bins" / original_fasta.name)
        manifest.append({"bin_id": original_id, "parent": parent, "k": "", "role": "parent", "child": ""})

        clusters: dict[str, set[str]] = {}
        for proposal in sorted(proposals, key=lambda row: int(row["k"])):
            assignment_path = Path(proposal["proposal_dir"]) / "candidate_children.tsv"
            proposal_clusters: dict[str, set[str]] = defaultdict(set)
            for assignment in read_tsv(assignment_path):
                proposal_clusters[assignment["clustername"]].add(assignment["contigname"])
            assigned = set().union(*proposal_clusters.values()) if proposal_clusters else set()
            if assigned != parent_ids:
                raise ValueError(
                    f"{parent}/k{proposal['k']}: candidate assignment differs from canonical parent "
                    f"(missing={len(parent_ids - assigned)}, extra={len(assigned - parent_ids)})"
                )
            for cluster, contigs in sorted(proposal_clusters.items()):
                if cluster in clusters:
                    raise ValueError(f"Duplicate child identifier: {cluster}")
                clusters[cluster] = contigs
                manifest.append({
                    "bin_id": cluster,
                    "parent": parent,
                    "k": int(proposal["k"]),
                    "role": "child",
                    "child": cluster.rsplit("__child", 1)[-1],
                })
        if clusters:
            with vamb.vambtools.Reader(parent_fasta) as reader:
                vamb.vambtools.write_bins(bins_dir, clusters, reader, maxbins=None)

    with (args.output / "bin_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["bin_id", "parent", "k", "role", "child"], delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest)
    if not pending:
        (args.output / "NO_STABLE_PROPOSALS").touch()
    print(f"materialized stable_proposals={len(pending)} validation_bins={len(manifest)} root={args.output}")


def build_marker_file(args: argparse.Namespace) -> None:
    parent_lines: dict[str, str] = {}
    with args.parent_lineage.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            bin_id, remainder = line.rstrip("\n").split("\t", 1)
            if not bin_id.endswith("__original"):
                raise ValueError(f"Unexpected parent lineage bin id: {bin_id}")
            parent_lines[bin_id.removesuffix("__original")] = remainder
    manifest = read_tsv(args.manifest)
    missing = sorted({row["parent"] for row in manifest} - parent_lines.keys())
    if missing:
        raise ValueError(f"Parents missing from CheckM lineage marker file: {missing}")
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite {args.output}")
    with args.output.open("w") as handle:
        handle.write("# [Lineage Marker File]\n")
        for row in manifest:
            handle.write(f"{row['bin_id']}\t{parent_lines[row['parent']]}\n")
    print(f"wrote parent-matched marker sets for {len(manifest)} bins to {args.output}")


def load_marker_hits(path: Path) -> dict[str, dict[str, list[tuple[str, str]]]]:
    result: dict[str, dict[str, list[tuple[str, str]]]] = {}
    with path.open() as handle:
        for line in handle:
            bin_id, encoded = line.rstrip("\n").split("\t", 1)
            genes = ast.literal_eval(encoded)
            markers: dict[str, list[tuple[str, str]]] = defaultdict(list)
            for gene_id, marker_hits in genes.items():
                contig = gene_id.rsplit("_", 1)[0]
                for marker_id in marker_hits:
                    markers[marker_id].append((gene_id, contig))
            result[bin_id] = dict(markers)
    return result


def load_qa(path: Path) -> dict[str, dict[str, str]]:
    rows = read_tsv(path)
    return {row["Bin Id"]: row for row in rows}


def report(args: argparse.Namespace) -> None:
    manifest = read_tsv(args.manifest)
    marker_hits = load_marker_hits(args.marker_stats)
    qa = load_qa(args.qa)
    children: dict[tuple[str, int], list[str]] = defaultdict(list)
    parent_ids: dict[str, str] = {}
    for row in manifest:
        if row["role"] == "parent":
            parent_ids[row["parent"]] = row["bin_id"]
        else:
            children[(row["parent"], int(row["k"]))].append(row["bin_id"])

    conflict_rows: list[dict[str, object]] = []
    proposal_rows: list[dict[str, object]] = []
    for (parent, k), child_ids in sorted(children.items()):
        parent_id = parent_ids[parent]
        parent_markers = marker_hits.get(parent_id, {})
        duplicate_markers = {marker: hits for marker, hits in parent_markers.items() if len(hits) > 1}
        resolved = unresolved = discordant = 0
        for marker, parent_locations in sorted(duplicate_markers.items()):
            child_locations = {child: marker_hits.get(child, {}).get(marker, []) for child in sorted(child_ids)}
            total = sum(len(hits) for hits in child_locations.values())
            maximum = max((len(hits) for hits in child_locations.values()), default=0)
            if total != len(parent_locations):
                status = "DISCORDANT_COPY_COUNT"
                discordant += 1
            elif maximum > 1:
                status = "UNRESOLVED_WITHIN_CHILD"
                unresolved += 1
            else:
                status = "RESOLVED_ACROSS_CHILDREN"
                resolved += 1
            conflict_rows.append({
                "parent": parent,
                "k": k,
                "marker_id": marker,
                "status": status,
                "parent_copies": len(parent_locations),
                "parent_genes": ";".join(gene for gene, _ in parent_locations),
                "parent_contigs": ";".join(contig for _, contig in parent_locations),
                "child_copies": ";".join(f"{child}:{len(hits)}" for child, hits in child_locations.items()),
                "child_genes": ";".join(
                    f"{child}:{','.join(gene for gene, _ in hits)}" for child, hits in child_locations.items() if hits
                ),
                "child_contigs": ";".join(
                    f"{child}:{','.join(contig for _, contig in hits)}" for child, hits in child_locations.items() if hits
                ),
            })
        parent_qa = qa[parent_id]
        child_qa = [qa[child] for child in child_ids]
        all_resolved = bool(duplicate_markers) and unresolved == 0 and discordant == 0
        decision = "PENDING_COMPLEMENTARY_QC" if all_resolved else "RETAIN_PARENT_MARKER_GATE"
        proposal_rows.append({
            "parent": parent,
            "k": k,
            "parent_completeness": parent_qa["Completeness"],
            "parent_contamination": parent_qa["Contamination"],
            "child_count": len(child_ids),
            "child_hq": sum(float(row["Completeness"]) >= 90 and float(row["Contamination"]) < 5 for row in child_qa),
            "child_usable": sum(float(row["Completeness"]) >= 50 and float(row["Contamination"]) < 10 for row in child_qa),
            "duplicated_parent_markers": len(duplicate_markers),
            "resolved_markers": resolved,
            "unresolved_markers": unresolved,
            "discordant_copy_counts": discordant,
            "decision": decision,
            "reason": (
                "Marker conflicts resolve under the parent marker set; requires child-quality and complementary GUNC/taxonomic QC"
                if all_resolved else
                "Parent duplicated-marker conflicts do not resolve cleanly under this split"
            ),
        })

    with args.conflicts.open("w", newline="") as handle:
        fields = ["parent", "k", "marker_id", "status", "parent_copies", "parent_genes", "parent_contigs", "child_copies", "child_genes", "child_contigs"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(conflict_rows)
    with args.output.open("w", newline="") as handle:
        fields = ["parent", "k", "parent_completeness", "parent_contamination", "child_count", "child_hq", "child_usable", "duplicated_parent_markers", "resolved_markers", "unresolved_markers", "discordant_copy_counts", "decision", "reason"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(proposal_rows)
    print(f"wrote proposals={len(proposal_rows)} conflicts={len(conflict_rows)}; no split was accepted")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="stage", required=True)
    materialize_parser = subparsers.add_parser("materialize")
    materialize_parser.add_argument("--decisions", type=Path, required=True)
    materialize_parser.add_argument("--parent-bins", type=Path, required=True)
    materialize_parser.add_argument("--output", type=Path, required=True)
    marker_parser = subparsers.add_parser("build-marker-file")
    marker_parser.add_argument("--parent-lineage", type=Path, required=True)
    marker_parser.add_argument("--manifest", type=Path, required=True)
    marker_parser.add_argument("--output", type=Path, required=True)
    report_parser = subparsers.add_parser("report")
    report_parser.add_argument("--manifest", type=Path, required=True)
    report_parser.add_argument("--marker-stats", type=Path, required=True)
    report_parser.add_argument("--qa", type=Path, required=True)
    report_parser.add_argument("--conflicts", type=Path, required=True)
    report_parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    {"materialize": materialize, "build-marker-file": build_marker_file, "report": report}[args.stage](args)


if __name__ == "__main__":
    main()
