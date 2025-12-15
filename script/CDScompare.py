# -*- coding: utf-8 -*-
"""
Compute the similarites between two or more structural annotations (GFF) of a same genome.
"""

import sys, os
from script.python_util.argparser import build_parser
from script.python_util.comparison import annotation_comparison
from script.python_util.multi import multicomp


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    gffs = args.gff_files

    if len(gffs) < 2:
        parser.error("At least two GFF files must be provided.")

    if args.pairing_mode == "best":
        print("\nPairing mode set to 'best'. Genes from both annotations within the same cluster "
              "are aligned using a pairwise alignment to find the best global gene pairing.\n")
    else:
        print("\nPairing mode set to 'all'. Within a cluster, all annotation1/annotation2 "
              "gene pairings are output as soon as their mRNA regions overlap.\n")

    mode_align = args.pairing_mode == "best"
    
    # Simple mode
    if len(gffs) == 2:
        ref_path, alt_path = gffs
        return annotation_comparison(ref_path, alt_path, args.out_dir, mode_align)

    # Multi mode
    ref_path = gffs[0]
    alt_paths = gffs[1:]
    return multicomp(ref_path, alt_paths, args.out_dir, mode_align)


if __name__ == "__main__":
    main()
