# -*- coding: utf-8 -*-

import csv
from pathlib import Path

from cdscompare.python_util.annotation import AnnotationPair, AnnotationSet


PAIRWISE_HEADER = (
    "seqid_strand",
    "cluster",
    "annot1_gene",
    "annot2_gene",
    "matches",
    "mismatches",
    "similarity_score",
    "annot1_start",
    "annot1_end",
    "annot2_start",
    "annot2_end",
    "annot1_mRNA",
    "annot2_mRNA",
    "C_NC_mismatch_zones",
    "RF_mismatch_zones",
    "C_NC_mismatches",
    "RF_mismatches",
    "annot1_mRNA_count",
    "annot2_mRNA_count",
)


def format_mismatch_zones(zones: list[int]) -> str:
    """Format mismatch coordinate pairs for CSV output."""
    return " ".join(
        f"[{zones[i]}//{zones[i + 1]}]"
        for i in range(0, len(zones), 2)
    )


def format_comparison_summary(
    stats: list[int],
    pair: AnnotationPair,
    seqid_strand: str | None = None,
) -> str:
    """Format pairwise comparison counts."""
    location = f" on {seqid_strand}" if seqid_strand else ""

    return (
        f"Comparison summary for {pair.ref.id} vs {pair.alt.id}{location}:\n"
        f"- reported gene pairs: {stats[0]}\n"
        f"- unpaired genes in annotation 1: {stats[1]}\n"
        f"- unpaired genes in annotation 2: {stats[2]}\n"
    )


def build_result_row(
    seqid_strand: str,
    result: dict,
) -> list[str | int]:
    """Build one detailed pairwise CSV row."""
    row: list[str | int] = [
        seqid_strand,
        result["cluster name"],
        result["reference"],
        result["alternative"],
    ]

    if not result["mismatch/match"]:
        return row + [
            "_",
            "_",
            f"{result['identity']:.2f}",
            result["reference start"],
            result["reference end"],
            result["alternative start"],
            result["alternative end"],
            result["reference mRNA"],
            result["alternative mRNA"],
            "_",
            "_",
            "_",
            "_",
            result["reference mRNA number"],
            result["alternative mRNA number"],
        ]

    matches, c_nc_mismatches, rf_mismatches = result["mismatch/match"]
    c_nc_zones, rf_zones = result["mismatch zones"]

    return row + [
        matches,
        c_nc_mismatches + rf_mismatches,
        f"{result['identity']:.2f}",
        result["reference start"],
        result["reference end"],
        result["alternative start"],
        result["alternative end"],
        result["reference mRNA"],
        result["alternative mRNA"],
        format_mismatch_zones(c_nc_zones),
        format_mismatch_zones(rf_zones),
        c_nc_mismatches,
        rf_mismatches,
        result["reference mRNA number"],
        result["alternative mRNA number"],
    ]


def update_comparison_stats(stats: list[int], result: dict) -> None:
    """Update pair and unpaired-gene counters from one result."""
    if result["mismatch/match"]:
        stats[0] += 1
        return

    if result["reference"] == "~":
        stats[2] += 1
        return

    stats[1] += 1


def write_results(
    all_results: dict,
    csv_path: Path,
    txt_path: Path,
    pair: AnnotationPair,
) -> None:
    """Write detailed pairwise results and comparison summaries."""
    full_stats = [0, 0, 0]

    csv_path.parent.mkdir(parents=True, exist_ok=True)

    with csv_path.open("w", newline="", encoding="utf-8") as results_file:
        csv_writer = csv.writer(results_file, lineterminator="\n")
        csv_writer.writerow(PAIRWISE_HEADER)

        for seqid_strand, clusters in all_results.items():
            seqid_stats = [0, 0, 0]

            for cluster in clusters:
                for result in cluster:
                    csv_writer.writerow(
                        build_result_row(seqid_strand, result)
                    )
                    update_comparison_stats(seqid_stats, result)

            summary = format_comparison_summary(
                seqid_stats,
                pair,
                seqid_strand,
            )
            print(f"\n{summary}", end="")

            full_stats = [
                total + current
                for total, current in zip(full_stats, seqid_stats)
            ]

    summary = format_comparison_summary(full_stats, pair)
    txt_path.write_text(summary, encoding="utf-8")
    print(f"\n{summary}", end="")


def write_multi_results(
    multi_results: list[dict],
    annotations: AnnotationSet,
    out_dir: Path,
) -> None:
    """Write the multi-comparison synthesis CSV file."""
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = annotations.synthesis_filename(out_dir)

    header = [f"{annotations.ref.id}_gene"]

    for annotation in annotations.alts:
        header.extend(
            [
                f"{annotation.id}_gene",
                f"{annotation.id}_similarity_score",
            ]
        )

    reference_genes: set[str] = set()

    for result in multi_results:
        reference_genes.update(result)

    with csv_path.open("w", newline="", encoding="utf-8") as results_file:
        csv_writer = csv.writer(results_file, lineterminator="\n")
        csv_writer.writerow(header)

        for reference_gene in sorted(reference_genes):
            row: list[str] = [reference_gene]

            for compared_result in multi_results:
                comparison = compared_result.get(reference_gene)

                if comparison is None:
                    row.extend(["~", "0.00"])
                    continue

                compared_gene, similarity = comparison
                row.extend([compared_gene, f"{similarity:.2f}"])

            csv_writer.writerow(row)
