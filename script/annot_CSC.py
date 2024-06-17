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
import intervals_utils as iu
import locus as lc
import read_files as rf
import pre_comparison as pc
import comparison as comp



## This function writes to a new 'results.csv' file the results of the 
# annotation comparison retrieved from the identities dictionary returned by 
# the annotation_match function
#
# @see annotation_match()
#
# @param results A list of list of dictionaries containing results of the 
# annotation comparison, as returned by annotation_match
#
# @param debug If True, triggers display of many messages intended for 
# debugging the program. Default is 'False'
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @remark Results are written in CSV ('Comma-Separated Values') format
def write_results(results, debug=False, verbose=False):
    
    # try to open the results dumping file
    try:
        results_file = open("./results/results.csv", "w")
    except FileNotFoundError:
        os.mkdir("./results/") # create 'results' subdirectory
        results_file = open("./results/results.csv", "w")
    
    results_file.write("Cluster name, Reference locus,Alternative locus,Comparison matches,Comparison mismatches,Identity score (%),Reference start, Reference end, Alternative start, Alternative end, Reference mRNA, Alternative mRNA, non-correspondance zones, reference mRNA number, alternative mRNA number\n")
        
    # annotation origin of each locus in the results
    # (first value : both,  second value : reference,  third value : alternative)
    locus_initial_annot = [0,0,0]    
    
    for cluster in results:
        for loc in cluster:
            # convert the mismatch zones so that commas don't modify 
            # the csv structure
            mismatch_zones = ""
            
            if loc["mismatch zones"][0] not in ["_", "?"]:
                if debug: print(f"mismatch zones = {loc['mismatch zones']}")
                for coords in loc["mismatch zones"]:
                    mismatch_zones += "[" + str(coords[0]) + "//" + str(coords[1]) + "] "
                            
            # if no comparison was done for the loci, write '~' instead of 
            # the comparison values
            if loc['mismatch/match'] == []:
                print(f"{loc['cluster name']}\t\t{loc['reference']}\t\t{loc['alternative']}\t\t\t_\t\t\t\t_")
                results_file.write(f"{loc['cluster name']},{loc['reference']},{loc['alternative']},_,_,{loc['identity']},{loc['reference start']},{loc['reference end']},{loc['alternative start']},{loc['alternative end']}, {loc['reference mRNA']}, {loc['alternative mRNA']}, _, {loc['reference mRNA number']}, {loc['alternative mRNA number']}\n")       
                if loc['reference'] == '~':
                    locus_initial_annot[2] += 1
                else:                       
                    locus_initial_annot[1] += 1
            else:
                print(f"{loc['cluster name']}\t\t{loc['reference']}\t\t{loc['alternative']}\t\t\t{loc['mismatch/match']}\t\t\t\t{loc['identity']}%")
                results_file.write(f"{loc['cluster name']},{loc['reference']},{loc['alternative']},{loc['mismatch/match'][1]},{loc['mismatch/match'][0]},{loc['identity']},{loc['reference start']},{loc['reference end']},{loc['alternative start']},{loc['alternative end']}, {loc['reference mRNA']}, {loc['alternative mRNA']}, {mismatch_zones}, {loc['reference mRNA number']}, {loc['alternative mRNA number']}\n")
                locus_initial_annot[0] += 2
                
    print(f"\nNumber of loci:\n- found in both annotations : {locus_initial_annot[0]}\n- found only in the reference : {locus_initial_annot[1]}\n- found only in the alternative : {locus_initial_annot[2]}\n")
    results_file.close()
    

## Main function of this program. Given a reference and alternative path, 
# gets the corresponding GFF files and compares the two annotations to return 
# their information about their loci's comparison
#
# @param ref_path Path of the GFF file describing the reference annotation
#
# @param alt_path Path of the GFF file describing the aternative annotation
#
# @param debug If True, triggers display of many messages intended for 
# debugging the program. Default is 'False'
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @param exon_mode Boolean indicating if the main comparison structures read
# from the file should be coding sequences (CDS, False) or exons (True).
# Default is 'False' (CDS comparison)
#
# @return Returns a list of lists of dictionaries describing the 
# comparison of the structure identity between the loci of each annotation 
def annotation_comparison(ref_path, alt_path, debug=False, verbose=False, create_strings=False, exon_mode=False):

    # get all annotation files and generate the annotation data structure
    ref_annotations = rf.get_gff_borders(ref_path, debug, verbose, exon_mode)
    alt_annotations = rf.get_gff_borders(alt_path, debug, verbose, exon_mode)
    
    # get the order of the loci of both annotations
    locus_order = pc.annotation_sort(ref_annotations, alt_annotations, debug, verbose)
    
    # construct clusters of overlapping loci
    cluster_list = pc.construct_clusters(ref_annotations, alt_annotations, locus_order, debug, verbose)
    
    results = []
    for cluster_id, cluster in cluster_list.items():
        results.append(comp.annotation_match(cluster, create_strings, debug, verbose))
        
    print("\nCluster name\tReference_Locus\t\tAlternative_Locus\t\tComparison[match/mismatch_EI/mismatch_RF]\t\tIdentity_Score\n")
    write_results(results, debug, verbose)
    
    return results
    

def usage():
    
    # displayed when '-h' or '--help' is given, or when an invalid script
    # call happens
    print("Syntax : path/to/main.py [ -h/--help -d/--debug -v/--verbose -o/--old_version ] [ -r/--reference <reference_file_path> ] [ -a/--alternative <alternative_file_path> ] ")
    

def main():
    
    # we retrieve all script call options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hdvoer:a:", ["help", "debug", "verbose", "old_version", "exon-mode", "reference=", "alternative="])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)
        
    # initialisation of the display parameters 
    #
    # debug: display lots of informations on internal function variables and 
    # mecanisms, not intended to be used by the end user
    #
    # verbose: display messages indicating which step the program is currently 
    # on, intended to be used when the program is called directly (not
    # integrated in a pipeline or workflow)
    debug = False
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
        elif o in ("-d", "--debug"):
            debug = True
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
    return annotation_comparison(ref_path, alt_path, debug, verbose, create_strings, exon_mode)
    
if __name__ == "__main__":
    main()
    #test_reverse()

    
def test_reverse():
    loc = lc.Locus(name="locus1", mRNAs={"chrblabla": [100, 200, 300, 400]}, start=100, end=400, direction="reverse")
    loc.reverse(600)
    print(loc.mRNAs)
    #FAIRE LA REINVERSION POUR LES ZONES DE MISMATCH

def test_get_reading_frame():
    cds_bounds_ref =[4,9,12,13]
    cds_bounds_alt=[1,6,9,10]
    cds_inter = [4,6,9,9]
    rf_ref=get_reading_frame(cds_bounds_ref, cds_inter, True, True)
    rf_alt=get_reading_frame(cds_bounds_alt, cds_inter, True, True)
    print(rf_ref)
    print(rf_alt)

def test2_get_reading_frame():
    cds_bounds_ref =[4,9,12,13]
    cds_bounds_alt=[4,13]
    cds_inter = [4,9,12,13]
    rf_ref=get_reading_frame(cds_bounds_ref, cds_inter, True, True)
    rf_alt=get_reading_frame(cds_bounds_alt, cds_inter, True, True)
    print(rf_ref)
    print(rf_alt)

def test3_get_reading_frame():
    cds_bounds_ref =[4,9,12,13]
    cds_bounds_alt=[5,9,12,13]
    cds_inter = [5,9,12,13]
    rf_ref=get_reading_frame(cds_bounds_ref, cds_inter, True, True)
    rf_alt=get_reading_frame(cds_bounds_alt, cds_inter, True, True)
    print(rf_ref)
    print(rf_alt)

def test_mrna_comp():
    cds_bounds_ref =[4,9,12,13]
    cds_bounds_alt=[1,6,9,10]
    intervals_ref = iu.OrderedIntervals(cds_bounds_ref, True);
    (matches, mismatches_EI, mismatches_RF, diff_EI, diff_RF) = compute_matches_mismatches_EI_RF (cds_bounds_ref, intervals_ref, cds_bounds_alt)
    print([matches, mismatches_EI, mismatches_RF])
    print( matches / (matches + mismatches_EI+mismatches_RF) * 100)

def test_mrna_comp2():
    cds_bounds_ref =[4,9,12,13]
    cds_bounds_alt=[4,13]
    intervals_ref = iu.OrderedIntervals(cds_bounds_ref, True);
    (matches, mismatches_EI, mismatches_RF, diff_EI, diff_RF) = compute_matches_mismatches_EI_RF (cds_bounds_ref, intervals_ref, cds_bounds_alt)
    print([matches, mismatches_EI, mismatches_RF])
    print( matches / (matches + mismatches_EI+mismatches_RF) * 100)


def test_mrna_comp3():
    cds_bounds_ref =[4,9,12,13]
    cds_bounds_alt=[5,9,12,13]
    intervals_ref = iu.OrderedIntervals(cds_bounds_ref, True);
    (matches, mismatches_EI, mismatches_RF, diff_EI, diff_RF) = compute_matches_mismatches_EI_RF (cds_bounds_ref, intervals_ref, cds_bounds_alt)
    print([matches, mismatches_EI, mismatches_RF])
    print( matches / (matches + mismatches_EI+mismatches_RF) * 100)

def test_mrna_comp4():
    cds_bounds_ref =[4,9,12,13]
    cds_bounds_alt=[5,9,11,13]
    intervals_ref = iu.OrderedIntervals(cds_bounds_ref, True);
    (matches, mismatches_EI, mismatches_RF, diff_EI, diff_RF) = compute_matches_mismatches_EI_RF (cds_bounds_ref, intervals_ref, cds_bounds_alt)
    print([matches, mismatches_EI, mismatches_RF])
    print( matches / (matches + mismatches_EI+mismatches_RF) * 100)
    