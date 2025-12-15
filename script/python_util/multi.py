#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from python_util.comparison import annotation_comparison
from python_util.io import write_multi_results


## Extracts the computed identity for each locus of the given result dictionary
# and add it to a new dictionary
#
# @param result Result dictionary, as returned by CDScompare's 
# annotation_compare() function
#
# @returns Returns a dictionary with the locus name as key, and the comparison
# identity as value
def result_to_dict(result):
    # key is the ref gene id, value is the tuple (alt_gene_id, identity)
    dict_result = {}
    for chrm_name, chrm in result.items():
        for cluster in chrm:
            for locus in cluster:
                ref_locus = locus['reference']
                alt_locus = locus['alternative']
                identity = locus['identity']
                if(ref_locus != "~"):
                    dict_result[ref_locus] = (alt_locus, identity)
    return dict_result




## Compares all alternative annotations given (alt_paths) to the reference 
# annotation (ref_path) by calling annotation_sort and appends the loci 
# identities to a list
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
# @returns Returns the list of all loci identities (dictionaries) of all 
# alternatives
def multicomp(ref_path, alt_paths, out_dir, mode_align):
    
    # list of all results dictionaries returned by CDScompare
    multi_results = []
    
    # for each alternative annotation given to the program, use CDScompare to 
    # compute results for the comparison with the reference, write the results
    # in a CSV file, and append the returned loci identities to a list
    for alt_path in alt_paths:
        multi_results.append(result_to_dict(annotation_comparison(ref_path, alt_path, out_dir, mode_align)))

        
    write_multi_results(multi_results, ref_path, alt_paths)
    
    return multi_results
    
