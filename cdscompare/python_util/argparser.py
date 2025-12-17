# -*- coding: utf-8 -*-

import argparse

def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Compare structural genome annotations (CDScompare).\n"
            "When several genes overlap, they are clustered together and compared pairwise.\n"
            "Results are written to a CSV file in the specified output directory.\n"
        ),
        epilog = (
            "Pairing modes:\n"
            "  best (default): Genes from both annotations within the same cluster are aligned using a pairwise alignment to find the best global gene pairing.\n"
            "  all: Within a cluster, all annotation1/annotation2 gene pairings are output as soon as their mRNA regions overlap.\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "gff_files",
        nargs="+",
        help="Two or more GFF files to compare."
    )

    parser.add_argument(
        "-d", "--out_dir",
        default="results",
        help="Output directory (default: results)"
    )

    parser.add_argument(
        "-p", "--pairing_mode",
        choices=["best", "all"],
        default="best",
        help="Pairing mode for CDS comparison."
    )

    return parser
