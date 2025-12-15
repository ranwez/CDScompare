#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 11:57:29 2024

@author: vetea, ranwez
"""

##@package CDScompare
# This script is used to compute the similarites between two or more structural
# annotations of a same genome.
# It expects as input the paths to the annotation files (in GFF format),
# displays the computed similarities between all annotation pairs, and creates
# a results CSV file detailing the loci comparisons between the annotations.

# import getopt
import sys
import os
from python_util.argparser import build_parser
from python_util.comparison import annotation_comparison
from python_util.multi import multicomp

script_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "python_util")
sys.path.append( script_dir )

def main():
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
