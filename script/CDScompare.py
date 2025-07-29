#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 11:57:29 2024

@author: vetea, ranwez
"""

##@package CDScompR
# This script is used to compute the similarites between two structural
# annotations of a same genome, one reference and one alternative annotation.
# It expects as input the paths to the annotation files (in GFF format),
# displays the computed similarities between all annotation pairs, and creates
# a results CSV file detailing the loci comparisons between the annotations

import getopt
import sys
import os

script_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "python_util")
sys.path.append( script_dir )

import comparison as comp
import cluster as cl
from locus import gff_to_cdsInfo, Locus, STRING_CACHE_REVERSE

def format_mismatch_zones(zones):
    """Formate une liste de coordonnées de mismatch en chaîne lisible"""
    return " ".join(f"[{zones[i]}//{zones[i+1]}]" for i in range(0, len(zones), 2))

def write_results(all_results, alt_name, out_dir):
    # gene stats [gene found in both; found only in ref; found only in alt]
    full_annotation_stat = [0, 0, 0]

    # try to open the results file (named after the alternative file name)
    filename=f"{out_dir}/{alt_name}.csv"
    try:
        results_file = open(filename, "w")
    except FileNotFoundError:
        os.mkdir(out_dir) # create 'results' subdirectory
        results_file = open(filename, "w")

    results_file.write("Chromosome,Cluster name,Reference locus,Alternative locus,Comparison matches,Comparison mismatches,Identity score (%),Reference start,Reference end,Alternative start,Alternative end,Reference mRNA,Alternative mRNA,Exon_intron (EI) non-correspondance zones,Reading frame (RF) non-correspondance zones,Exon_Intron (EI) mismatches,Reading Frame (RF) mismatches,reference mRNA number,alternative mRNA number\n")

    for dna_mol, results in all_results.items():
        chromosome_strand_stat = [0,0,0]

        for cluster in results:
            for loc in cluster:
                if loc['mismatch zones'] not in ["_", "?"]:
                    # convert mismatch zones so commas don't modify the CSV output
                    mismatch_EI = format_mismatch_zones(loc['mismatch zones'][0])
                    mismatch_RF = format_mismatch_zones(loc['mismatch zones'][1])

                else:
                    mismatch_EI = loc['mismatch zones']
                    mismatch_RF = loc['mismatch zones']

                if loc['mismatch/match'] == []:
                    results_file.write(f"{dna_mol},{loc['cluster name']},{loc['reference']},{loc['alternative']},_,_,{loc['identity']},{loc['reference start']},{loc['reference end']},{loc['alternative start']},{loc['alternative end']},{loc['reference mRNA']},{loc['alternative mRNA']},_,_,_,_,{loc['reference mRNA number']},{loc['alternative mRNA number']}\n")
                    if loc['reference'] == '~':
                        chromosome_strand_stat[2] += 1
                    else:
                        chromosome_strand_stat[1] += 1
                else:
                    results_file.write(f"{dna_mol},{loc['cluster name']},{loc['reference']},{loc['alternative']},{loc['mismatch/match'][0]},{loc['mismatch/match'][1]+loc['mismatch/match'][2]},{loc['identity']},{loc['reference start']},{loc['reference end']},{loc['alternative start']},{loc['alternative end']},{loc['reference mRNA']},{loc['alternative mRNA']},{mismatch_EI},{mismatch_RF},{loc['mismatch/match'][1]},{loc['mismatch/match'][2]},{loc['reference mRNA number']},{loc['alternative mRNA number']}\n")
                    chromosome_strand_stat[0] += 1

        print(f"\nNumber of loci of {dna_mol}:\n- found in both annotations : {chromosome_strand_stat[0]}\n- found only in the reference : {chromosome_strand_stat[1]}\n- found only in the alternative : {chromosome_strand_stat[2]}\n")
        full_annotation_stat = [sum(x) for x in zip(full_annotation_stat, chromosome_strand_stat)]

    results_file.close()

    results_file = open(f"{out_dir}/{alt_name}.txt", "w")
    results_file.write(f"\nNumber of loci (whole data):\n- found in both annotations : {full_annotation_stat[0]}\n- found only in the reference : {full_annotation_stat[1]}\n- found only in the alternative : {full_annotation_stat[2]}\n")
    results_file.close()
    print(f"\nNumber of loci (whole data):\n- found in both annotations : {full_annotation_stat[0]}\n- found only in the reference : {full_annotation_stat[1]}\n- found only in the alternative : {full_annotation_stat[2]}\n")

def build_cluster_list_from_Locus(ref_loci, alt_loci, dna_mol):
    """
    Build clusters of loci from reference and alternative annotations.
    
    Args:
        read_ref: Dictionary mapping chromosome_strand to list of Locus objects from reference annotation
        read_alt: Dictionary mapping chromosome_strand to list of Locus objects from alternative annotation
        dna_mol: Chromosome_strand key to process
        
    Returns:
        List of Cluster objects
    """
    cluster_list = []
    
    # Get the max end positions of the loci from both annotations
    sentinel = 1 + max(
        max(locus.end for locus in ref_loci) if ref_loci else 0,
        max(locus.end for locus in alt_loci) if alt_loci else 0
    )

    ref_sentinel = Locus.sentinel("sentinel_ref", sentinel)
    alt_sentinel = Locus.sentinel("sentinel_alt", sentinel)

    ref_loci.append(ref_sentinel)
    alt_loci.append(alt_sentinel)
    
    ref_i, alt_i = 0, 0
    ref_locus = ref_loci[ref_i]
    alt_locus = alt_loci[alt_i]
    
    cluster_max = -1
    cluster_loci = {"ref": [], "alt": []}
    cluster_id = 0
    
    while ref_locus.start < sentinel or alt_locus.start < sentinel:
        if ref_locus.start <= alt_locus.start:
            locus = ref_locus
            type = "ref"
            ref_locus, ref_i = ref_loci[ref_i+1], ref_i + 1
        else:
            locus = alt_locus
            type = "alt"
            alt_locus, alt_i = alt_loci[alt_i+1], alt_i + 1
        
        if locus.start > cluster_max:
            if cluster_max != -1:
                cluster_name = f"{dna_mol}_{cluster_id}"
                cluster_list.append(cl.Cluster(cluster_name, cluster_loci["ref"], cluster_loci["alt"], cluster_max))
            cluster_id += 1
            cluster_loci = {"ref": [], "alt": []}
            cluster_max = -1
        
        cluster_max = max(cluster_max, locus.end)
        cluster_loci[type].append(locus)
    
    if cluster_max != -1:
        cluster_name = f"{dna_mol}_{cluster_id}"
        cluster_list.append(cl.Cluster(cluster_name, cluster_loci["ref"], cluster_loci["alt"]  , cluster_max))

    ref_loci.clear()
    alt_loci.clear()
    return cluster_list

def annotation_comparison(ref_path:str, alt_path:str, out_dir:str, mode_align:bool):
    """Compare two GFF annotations and write results to a file."""
    read_ref= gff_to_cdsInfo(ref_path)
    read_alt= gff_to_cdsInfo(alt_path)
    all_results = {}
    reverse_str="_"+ STRING_CACHE_REVERSE
    dna_mols = list(read_ref.keys() | read_alt.keys())
    dna_mols.sort()
    for dna_mol in dna_mols:
        clusters = build_cluster_list_from_Locus(read_ref[dna_mol], read_alt[dna_mol], dna_mol) 
        results = [None] * len(clusters)
        for i, cluster in enumerate(clusters):
            results[i]=comp.annotation_match(cluster, dna_mol.endswith(reverse_str), mode_align)
        all_results[dna_mol] = results

    alt_name = (os.path.basename(alt_path)).split(".")[0]
    write_results(all_results, alt_name, out_dir)

    return all_results


def usage():
    print("Syntax : path/to/CDScompare.py [ -h/--help] [ -r/--reference <reference_file_path> ] [ -a/--alternative <alternative_file_path> ] [-d/--out_dir <output_directory>] [-p/--pairing_mode <best/all>] ")
    print("When several genes overlap, they are clustered together and compared pairwise. The results are written to a CSV file in the specified output directory.\n")
    print("In pairing mode 'best' (default), genes of the reference and alternative annotations within the same cluster are aligned using a pairwise alignment, to find the best global gene pairing. \n")
    print("In pairing mode 'all', all gene pairings, of the reference and alternative annotations within the same cluster, will be outputted as soon as their mRNA regions overlap.\n")

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hvoer:a:d:p:", ["help", "reference=", "alternative=", "out_dir=", "pairing_mode="])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)

    out_dir = "results"
    mode_align=True
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-r", "--reference"):
            ref_path = a
        elif o in ("-a", "--alternative"):
            alt_path = a
        elif o in ("-d", "--out_dir"):
            out_dir = a
        elif o in ("-p", "--pairing_mode"):
            if a not in ["best", "all"]:
                print("Invalid pairing mode. Use 'best' or 'all'.")
                sys.exit(2)
            pairing_mode = a
            if a == "all":
                mode_align = False
            if pairing_mode == "best":
                print("Pairing mode set to 'best'. Whithin cluster of overlapping genes, a pairwise alignment will be used to identify optimal gene pairings.")
            elif pairing_mode == "all":
                print("Pairing mode set to 'all'. Whithin cluster of overlapping genes, all ref/alt gene pairings will be outputted as soon as their mRNA regions overlap.")
        else:
            assert False, "unhandled option"

    return annotation_comparison(ref_path, alt_path, out_dir, mode_align)

if __name__ == "__main__":
    main()
