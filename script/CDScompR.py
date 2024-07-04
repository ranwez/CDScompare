#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 11:57:29 2024

@author: vetea
"""

##@package main
# This script is used to compute the distance between two structural 
# annotations of a same genome, one reference and one alternative annotation. 
# It expects as input the paths to the annotation files (in GFF format), 
# displays the computed distances between all annotation pairs, and returns a 
# dictionary of lists of lists detailing the matchs/mismatchs between the 
# two annotations' structure string


import getopt
import sys
import os
import read_files as rf
import pre_comparison as pc
import comparison as comp
import time

## This function writes to a new 'results.csv' file the results of the 
# annotation comparison retrieved from the identities dictionary returned by 
# the annotation_match function
#
# @see annotation_match()
#
# @param results A dictionary of list of list of dictionaries containing 
# results of the annotation comparison, as returned by annotation_match
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @see annotation_match()
#
# @remark Results are written in CSV ('Comma-Separated Values') format
def write_results(all_results, verbose=False):
    # annotation origin of each locus in the results for all chromosomes
    # (first value : both,  second value : reference,  third value : alternative)
    final_locus_annot = [0, 0, 0]
    
    # try to open the results file
    try:
        results_file = open("./results/results.csv", "w")
    except FileNotFoundError:
        os.mkdir("./results/") # create 'results' subdirectory
        results_file = open("./results/results.csv", "w")
        
    results_file.write("Chromosome, Cluster name, Reference locus,Alternative locus,Comparison matches,Comparison mismatches,Identity score (%),Reference start, Reference end, Alternative start, Alternative end, Reference mRNA, Alternative mRNA, Exon_intron (EI) non-correspondance zones, Reading frame (RF) non-correspondance zones, Exon_Intron (EI) mismatches, Reading Frame (RF) mismatches, reference mRNA number, alternative mRNA number\n")
    
    for dna_mol, results in all_results.items():
        print(f"\n**************** Results for chromosome {dna_mol} ****************\n")
            
        # annotation origin of each locus in the results
        # (first value: both,  second: reference,  third: alternative)
        locus_initial_annot = [0,0,0]
        
        for cluster in results:
            for loc in cluster:
                if loc['mismatch zones'] not in ["_", "?"]:
                    # convert mismatch zones so commas don't modify the CSV output
                    mismatch_EI = ""
                    mismatch_RF = ""
                    for i in range(0, len(loc['mismatch zones'][0]), 2):
                        mismatch_EI += "[" + str(loc['mismatch zones'][0][i]) + "//" + str(loc['mismatch zones'][0][i+1]) + "] "
                    for i in range(0, len(loc['mismatch zones'][1]), 2):
                        mismatch_RF += "[" + str(loc['mismatch zones'][1][i]) + "//" + str(loc['mismatch zones'][1][i+1]) + "] "
                    
                else:
                    mismatch_EI = loc['mismatch zones']
                    mismatch_RF = loc['mismatch zones']
                
                # if no comparison was done for the loci, write '~' instead of 
                # the comparison values
                if loc['mismatch/match'] == []:
                    print(f"{loc['cluster name']}\t\t{loc['reference']}\t\t{loc['alternative']}\t\t\t_\t\t\t\t_")
                    results_file.write(f"{dna_mol}, {loc['cluster name']},{loc['reference']},{loc['alternative']},_,_,{loc['identity']},{loc['reference start']},{loc['reference end']},{loc['alternative start']},{loc['alternative end']}, {loc['reference mRNA']}, {loc['alternative mRNA']}, _, _, _, _, {loc['reference mRNA number']}, {loc['alternative mRNA number']}\n")       
                    if loc['reference'] == '~':
                        locus_initial_annot[2] += 1
                    else:                       
                        locus_initial_annot[1] += 1
                else:
                    print(f"{loc['cluster name']}\t\t{loc['reference']}\t\t{loc['alternative']}\t\t\t{loc['mismatch/match']}\t\t\t\t{loc['identity']}%")
                    results_file.write(f"{dna_mol}, {loc['cluster name']},{loc['reference']},{loc['alternative']},{loc['mismatch/match'][0]},{loc['mismatch/match'][1]+loc['mismatch/match'][2]},{loc['identity']},{loc['reference start']},{loc['reference end']},{loc['alternative start']},{loc['alternative end']}, {loc['reference mRNA']}, {loc['alternative mRNA']}, {mismatch_EI}, {mismatch_RF}, {loc['mismatch/match'][1]}, {loc['mismatch/match'][2]}, {loc['reference mRNA number']}, {loc['alternative mRNA number']}\n")
                    locus_initial_annot[0] += 2
                    
        print(f"\nNumber of loci of chromosome {dna_mol}:\n- found in both annotations : {locus_initial_annot[0]}\n- found only in the reference : {locus_initial_annot[1]}\n- found only in the alternative : {locus_initial_annot[2]}\n")
        # add locus origin counts of the chromosome to final counts for all
        final_locus_annot = [sum(x) for x in zip(final_locus_annot, locus_initial_annot)]
        
    results_file.close()
    print(f"\nNumber of loci (whole data):\n- found in both annotations : {final_locus_annot[0]}\n- found only in the reference : {final_locus_annot[1]}\n- found only in the alternative : {final_locus_annot[2]}\n")
    

## Main function of this program. Given a reference and alternative path, 
# gets the corresponding GFF files and compares the two annotations to return 
# their information about their loci's comparison
#
# @param ref_path Path of the GFF file describing the reference annotation
#
# @param alt_path Path of the GFF file describing the aternative annotation
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @param create_strings Boolean indicating wether to use the 'old' comparison
# function (old_compare_loci, 'True') or the new one (compare_loci, 'False')
#
# @param exon_mode Boolean indicating if the main comparison structures read
# from the file should be coding sequences (CDS, False) or exons (True).
# Default is 'False' (CDS comparison)
#
# @see compare_loci()
#
# @see old_compare_loci()
#
# @return Returns a list of lists of dictionaries describing the 
# comparison of the structure identity between the loci of each annotation
def annotation_comparison(ref_path, alt_path, verbose=False, create_strings=False, exon_mode=False):

    # get all annotation files and generate the annotation data structure
    ref_annotations = rf.get_gff_borders(ref_path, verbose, exon_mode)
    alt_annotations = rf.get_gff_borders(alt_path, verbose, exon_mode)

    # get the order of the loci of both annotations
    all_locus_order = {}
    for dna_mol in ref_annotations.keys():
        if verbose:
            print(f"Constructing the locus order list of the chromosome {dna_mol}")
        locus_order = pc.annotation_sort(ref_annotations[dna_mol], alt_annotations[dna_mol], verbose)
        all_locus_order[dna_mol] = locus_order
    
    # construct clusters of overlapping loci
    all_cluster_list = {}
    for dna_mol in all_locus_order.keys():
        cluster_list = pc.construct_clusters(ref_annotations[dna_mol], alt_annotations[dna_mol], all_locus_order[dna_mol], verbose)
        all_cluster_list[dna_mol] = cluster_list
    
    all_results = {}
    for dna_mol in all_cluster_list.keys():
        results = []
        for cluster_id, cluster in all_cluster_list[dna_mol].items():
            results.append(comp.annotation_match(cluster, create_strings, verbose))
        all_results[dna_mol] = results
        
    print("\nCluster name\tReference_Locus\t\tAlternative_Locus\t\tComparison[match/mismatch_EI/mismatch_RF]\t\tIdentity_Score\n")
    write_results(all_results, verbose)
    
    return all_results
    

def usage():
    
    # displayed when '-h' or '--help' is given, or when an invalid script
    # call happens
    print("Syntax : path/to/main.py [ -h/--help -v/--verbose -o/--old_version ] [ -r/--reference <reference_file_path> ] [ -a/--alternative <alternative_file_path> ] ")
    

def main():
    
    # we retrieve all script call options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hvoer:a:", ["help", "verbose", "old_version", "exon-mode", "reference=", "alternative="])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)
        
    # initialisation of the display parameters 
    #
    # verbose: display messages indicating which step the program is currently 
    # on, intended to be used when the program is called directly (not
    # integrated in a pipeline or workflow)
    verbose = False
    
    # boolean indicating which version of the file reading function to use: 
    # False = read coding sequences (CDS) from the given files, True = read
    # the exons from the given files
    exon_mode = False  
    
    # boolean indicating which version of the program to use: False = new 
    # version without any structure string creation, True = 'old' version with 
    # creation of structure strings to compare the loci of the annotations
    create_strings = False
    
    # we retrieve the values given for each parameter
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-v", "--verbose"):
            verbose = True
        elif o in ("-o", "--old_version"):
            create_strings = True
        elif o in ("-e", "--exon_mode"):
            exon_mode = True
        elif o in ("-r", "--reference"):
            ref_path = a
        elif o in ("-a", "--alternative"):
            alt_path = a
        else:
            assert False, "unhandled option"
            
    # call of the annotation_comparison function
    return annotation_comparison(ref_path, alt_path, verbose, create_strings, exon_mode)
    
if __name__ == "__main__":
    #annotation_comparison("../data/real_data/annot_best.gff", "../data/real_data/TRITD_clean.gff3", False, False, False)
    main()
    
    