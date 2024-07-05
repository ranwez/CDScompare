#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 13:18:44 2024

@author: vetea
"""

##@package CDSmulticompR
# This script is used to compute the distance between multiple structural 
# annotations of a same genome, one reference and one or more alternative 
# annotations. It expects as input the paths to the annotation files (in GFF 
# format), displays the computed distances between all annotation pairs, and 
# creates results CSV files for each alternative and one global synthesis file

import CDScompR as cc
import getopt
import sys
import os


## Writes the results returned by the function multicomp into a results 
# synthesis CSV file detailing the loci identity for each alternative
#
# @see multicomp()
#
# @param multi_results List of results dictionaries, as returned by multicomp
#
# @param ref_path Path to the reference annotation GFF file
def write_multi_results(multi_results, ref_path):
    
    ref_name = os.path.basename(ref_path).split(".")[0]
    
    # get each reference key of the first dictionary (since the reference is 
    # the same for all comparisons) and use it to retrieve the corresponding
    # comparisons for all result dictionaries, then write them in a synthesis
    # CSV file
    
    # try to open the results file
    try:
        results_file = open(f"./results/synthesis_{ref_name}.csv", "w")
    except FileNotFoundError:
        os.mkdir("./results/") # create 'results' subdirectory
        results_file = open(f"./results/synthesis_{ref_name}.csv", "w")
            
    # write the header    
    
    header = "Reference_locus"
    
    for alt in multi_results:
        alt_name = alt[1]
        header += f",{alt_name} locus,{alt_name} identity"
        
    results_file.write(header+"\n")
    
    # write the results
    
    for chrm_name, chrm in multi_results[0][0].items():
        cluster_id = 0
        for cluster in chrm:
            locus_id = 0
            for locus in cluster:
                ref_locus = locus['reference']
                line = f"{ref_locus}"
                
                for alt in multi_results:
                    # get name and identity of the associated alternative locus
                    alt_locus = alt[0][chrm_name][cluster_id][locus_id]['alternative']
                    alt_ident = alt[0][chrm_name][cluster_id][locus_id]['identity']
                    line += f",{alt_locus},{alt_ident}"
                locus_id += 1
                results_file.write(line+"\n")
                
            cluster_id += 1
                    
                    
    results_file.close()


## Compares all alternative annotations given (alt_paths) to the reference 
# annotation (ref_path) by calling annotation_sort and appends the results 
# dictionaries to a list
#
# @see annotation_sort()
#
# @param ref_path Path to the reference annotation file
#
# @param alt_path Path to the alternative annotation file
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
# @returns Returns the list of all results dictionaries of all alternatives
def multicomp(ref_path, alt_paths, verbose, create_strings, exon_mode):
    
    # list of all results dictionaries returned by CDScompR
    multi_results = []
    
    # for each alternative annotation given to the program, use CDScompR to 
    # compute results for the comparison with the reference, write the results
    # in a CSV file, and append the returned results dictionary to a list
    for alt in alt_paths:
        alt_name = os.path.basename(alt).split(".")[0]
        multi_results.append((cc.annotation_comparison(ref_path, alt, verbose, create_strings, exon_mode), alt_name))
        
    write_multi_results(multi_results, ref_path)
    
    return multi_results
    

def usage():
    
    # displayed when '-h' or '--help' is given, or when an invalid script
    # call happens
    print("Syntax : path/to/CDSmulticompR.py [ -h/--help -v/--verbose -o/--old_version ] [ -r/--reference <reference_file_path> ] [ -a/--alternative <alternative_file_path> -a/--alternative <other_alternative_file_path> -a/--alternative <...> ] ")
    

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
    
    # List of file paths to the alternative annotations files to compare
    # to the reference
    alt_paths = []
    
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
            alt_paths.append(a)
        else:
            assert False, "unhandled option"
            
            
    # call to the function responsible for the comparison of all files
    multicomp(ref_path, alt_paths, verbose, create_strings, exon_mode)
        
    
if __name__ == "__main__":
    main()