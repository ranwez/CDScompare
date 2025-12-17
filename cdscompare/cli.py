# -*- coding: utf-8 -*-
"""
Compute the similarites between two or more structural annotations (GFF) of a same genome.
"""
from pathlib import Path
from cdscompare.python_util.argparser import build_parser
from cdscompare.python_util.comparison import annotation_comparison
from cdscompare.python_util.multi import multicomp
from cdscompare.python_util.annotation import AnnotationSet

def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    gffs = [Path(p) for p in args.gff_files]
    out_dir = Path(args.out_dir)

    if len(gffs) < 2:
        parser.error("At least two GFF files must be provided.")

    if args.pairing_mode == "best":
        print("\nPairing mode set to 'best'. Genes from both annotations within the same cluster "
              "are aligned using a pairwise alignment to find the best global gene pairing.\n")
    else:
        print("\nPairing mode set to 'all'. Within a cluster, all annotation1/annotation2 "
              "gene pairings are output as soon as their mRNA regions overlap.\n")

    mode_align = args.pairing_mode == "best"
    
    annotations = AnnotationSet.from_paths(
        ref_path=gffs[0],
        alt_paths=gffs[1:],
    )

    # Simple mode
    if len(annotations.alts) == 1:
        pair = annotations.pairs()[0]
        annotation_comparison(pair, out_dir, mode_align)
        return

    # Multi mode
    multicomp(annotations, out_dir, mode_align)


if __name__ == "__main__":
    main()
