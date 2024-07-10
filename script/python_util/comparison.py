#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 09:38:04 2024

@author: vetea
"""

import intervals_utils as iu
import pre_comparison as pc


## reverses the coordinates of the mismatch zones returned by the functions
# compare_loci
#
# @see compare_loci()
#
# @param mismatch_zones Comparison mismatch zones as returned by compare_loci 
# (tuple of lists)
#
# @param cluster_end End coordinate of the last locus of the cluster
#
# @returns Returns a tuple of lists corresponding to the reversed zones
def reverse_coord(mismatch_zones, cluster_end):
    new_list_EI = []
    new_list_RF = []
    
    if mismatch_zones[0] != []:
        for bound in reversed(mismatch_zones[0]):
            new_list_EI.append(cluster_end-bound)
        
    if mismatch_zones[1] != []:
        for bound in reversed(mismatch_zones[1]):
            new_list_RF.append(cluster_end-bound)
            
    return (new_list_EI, new_list_RF)

## This function compares two annotations' loci returned by the function 
# get_gff_borders and creates for each pair of reference-alternative mRNAs 
# a comparison list detailing the identities and differences between the two
# annotations's codon position structure. It returns a tuple containing the 
# comparison list giving the highest identity, the computed identity level, 
# the list of mismatch areas indentified by the comparison, and the mRNAs
# with the highest identity value
#
# @see get_gff_borders()
#
# @param ref_locus The reference annotation's locus (Locus class instance)
#
# @param alt_locus The alternative annotation's locus (Locus class instance)
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @returns Returns a tuple containing a list indicating the number of codon 
# position mismatches (first position) and matches (second position) between 
# the two border lists giving the maximum indentity between all mRNAs, the 
# identity level computed from this list, the comaprison areas
# producing mismatches, and the compared mRNA names
#
# @remark This function doesn't expect any annotation to be a 'reference'
def compare_loci(ref_locus, alt_locus, verbose=False):
    final_comparison = [0,0,0]
    final_identity = 0.0
    final_mismatch_zones = []
    
    # if the loci are on different strands, return 'null' values
    if ref_locus.direction != alt_locus.direction:
        return ('_', 0.0, '_', '_', '_')
    
    # if the loci don't overlap, return 'null' values
    overlap = (alt_locus.start <= ref_locus.end <= alt_locus.end) or (ref_locus.start <= alt_locus.end <= ref_locus.end)
    if not overlap:
        return('_', 0.0, '_', '_', '_')
    
    for mRNA_ref_id, mRNA_ref in ref_locus.mRNAs.items():
        intervals_ref = iu.OrderedIntervals(mRNA_ref, True);
    
        for mRNA_alt_id, mRNA_alt in alt_locus.mRNAs.items():
            (matches, mismatches_EI, mismatches_RF, diff_EI, diff_RF) = compute_matches_mismatches_EI_RF(mRNA_ref, intervals_ref, mRNA_alt, verbose)
            identity = matches / (matches + mismatches_EI + mismatches_RF)
        
            # for each mRNA, we test wether the computed identity is higher 
            # than for the preceding mRNAs, to retrieve the highest identity
            if identity > final_identity:
                final_comparison = [matches, mismatches_EI, mismatches_RF]
                final_identity = identity
                final_mismatch_zones = (diff_EI, diff_RF)
                final_ref_mRNA = mRNA_ref_id
                final_alt_mRNA = mRNA_alt_id 
            # if all mRNAs comparisons return 0% identity, we still want 
            # mismatch values to be returned
            elif identity == 0.0 and mismatches_EI + mismatches_RF > final_comparison[1] + final_comparison[2]:
                final_comparison = [matches, mismatches_EI, mismatches_RF]
                final_mismatch_zones = (diff_EI, diff_RF)
                final_ref_mRNA = mRNA_ref_id
                final_alt_mRNA = mRNA_alt_id 
    
    # return the highest mRNA identity between the locus of each annotation
    final_identity = round(final_identity * 100, 1)
    return (final_comparison, final_identity, final_mismatch_zones, final_ref_mRNA, final_alt_mRNA)


## Computes the number of locus positions matching (both in CDS and same codon 
# position) or mismatching (one not in CDS or different codon position) in the
# given intervals
#
# @param mRNA_ref CDS position intervals of the reference (as a list)
#
# @param intervals_ref Instance of class 'OrderedIntervals' describing the 
# position intervals of the reference
#
# @param mRNA_alt CDS position intervals of the reference (as a list)
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @see OrderedIntervals
#
# @returns Returns a tuple of the number of positions matching, mismatching
# for 'EI', mismatching for 'RF', the 'EI' mismatch zones, and the 'RF' 
# mismatch zones
#
# @remark 'EI' stands for 'Exon-Intron' and indicates one of the annotations
# is not in a CDS for the current position; 'RF' stands for 'Reading Frame' and
# indicates the two annotations are not in the same codon position at the 
# current position
#
# @remark This function was written by Vincent Ranwez
def compute_matches_mismatches_EI_RF(mRNA_ref, intervals_ref, mRNA_alt, verbose):   
    intervals_alt = iu.OrderedIntervals(mRNA_alt, True);
    
    # get intersection and union of both intervals lists
    inter_mrna = intervals_ref.intersection(intervals_alt);
    union_mrna = intervals_ref.union(intervals_alt);
    
    # get exon and intron (EI) mismatches and mismatches zones (simple 
    # symmetric difference)
    mismatches_EI=union_mrna.total_length()-inter_mrna.total_length();
    diff_EI=intervals_ref.symmetric_difference(intervals_alt);
    
    matches=0 
    mismatches_RF=0; # reading frame (RF) mismatches
    diff_RF=[]; # reading frame (RF) mismatches zones
    
    # get intervals of intersection with their upper bounds
    inter_mrna_bounds=inter_mrna.get_intervals_with_included_ub();
    
    # get reading frames for each interval of the intersection
    rf_ref=pc.get_reading_frame(mRNA_ref, inter_mrna_bounds, True)
    rf_alt=pc.get_reading_frame(mRNA_alt, inter_mrna_bounds, True)
    
    # for each intersection interval of the reference and alternative, compare 
    # the reading frames to determine match or RF mismatch
    interval_id=0;
    for interval_id in range(0, len(rf_ref)):
        interval_lg= inter_mrna_bounds[2*interval_id+1]-inter_mrna_bounds[2*interval_id]+1;
        if rf_ref[interval_id] != rf_alt[interval_id]:
            mismatches_RF+=interval_lg
            diff_RF.append(inter_mrna_bounds[2*interval_id])
            diff_RF.append(inter_mrna_bounds[2*interval_id+1]);
        else:
            matches+=interval_lg
    return (matches, mismatches_EI, mismatches_RF, diff_EI.intervals, diff_RF)


## This function expects two loci corresponding to two annotations of the same 
# genome, creates and compares their structure strings (create_vectors function) 
#
# @see create_vectors()
#
# @param borders_vector_ref Instance of class 'Locus' of the reference locus
#
# @param borders_vector_alt Instance of class 'Locus' of the alternative locus
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @returns Returns a tuple containing the list of mismatches (first value) 
# and matches (second value) the computed string identity, and the
# compared mRNA names
#
# @remark This function doesn't expect any annotation to be a 'reference'
def old_compare_loci(ref_locus, alt_locus, verbose=False):
    if verbose:
        print(f"\n\n**************** comparing loci {ref_locus.name} of reference and {alt_locus.name} of alternative ****************")
    final_comparison = [0,0,0]
    final_identity = 0.0
    
    # if the loci are on different strands, return 'null' values
    if ref_locus.direction != alt_locus.direction:
        return ('_', 0.0, "_", "_")
    
    # if the loci don't overlap, return 'null' values
    overlap = (alt_locus.start <= ref_locus.end <= alt_locus.end) or (ref_locus.start <= alt_locus.end <= ref_locus.end)
    if not overlap:
        return('_', 0.0, "_", "_")
    
    for mRNA_ref_id, mRNA_ref in ref_locus.mRNAs.items():
        for mRNA_alt_id, mRNA_alt in alt_locus.mRNAs.items():
            if verbose:
                print(f"**************** comparing mRNA {mRNA_ref_id} of reference locus {ref_locus.name} and mRNA {mRNA_alt_id} of alternative locus {alt_locus.name} ****************")
            
            # we create the structure string for each mRNA
            start_ref, vector_ref = pc.create_vectors(mRNA_ref, verbose)
            start_alt, vector_alt = pc.create_vectors(mRNA_alt, verbose)
            comparison = [0,0,0]
                
            # we assign to 'minimum string' and 'maximum string' the reference
            # and alternative strings
            if start_ref <= start_alt:
                vector_min = vector_ref
                vector_max = vector_alt
                start_min = start_ref
                start_max = start_alt
            else:
                vector_min = vector_alt
                vector_max = vector_ref
                start_min = start_alt
                start_max = start_ref
            diff = start_max - start_min
            
            # for each position in the string of maximum length + position diff
            for i in range(max(len(vector_ref), len(vector_alt)) + diff):
                    
                # identitify the presence of a CDS at each position
                try:
                    min_in_cds = vector_min[i] in ["1", "2", "3"]    
                except IndexError:
                    min_in_cds = False
                
                # to prevent looping in case of negative 'i-diff', assign 
                # 'False' to max_in_cds if 'i-diff' is negative
                if i-diff >= 0:
                    try:
                        max_in_cds = vector_max[i-diff] in ["1", "2", "3"]
                    except IndexError:
                        max_in_cds = False
                else:
                    max_in_cds = False
                    
                # if both are in a cds and have the same codon position,
                # identify as a 'match'
                if min_in_cds and max_in_cds and vector_min[i] == vector_max[i-diff]:
                    comparison[0] += 1
                    
                # if both are in a cds but don't have the same codon
                # position, identify as a 'mismatch'
                elif min_in_cds and max_in_cds and vector_min[i] != vector_max[i-diff]:
                    comparison[2] += 1
                    
                # if only one is in a cds, identify as a 'mismatch'
                elif min_in_cds or max_in_cds:
                    comparison[1] += 1
                    
                # the cas of both not being in a cds is ignored
            identity = comparison[0] / (comparison[0] + comparison[1] + comparison[2])
            
            # for each mRNA, we test wether the computed identity is higher 
            # than for the preceding mRNAs, to retrieve the highest identity
            if identity > final_identity:
                final_comparison = comparison
                final_identity = identity
                final_ref_mRNA = mRNA_ref_id
                final_alt_mRNA = mRNA_alt_id 
            # if all mRNAs comparisons return 0% identity, we still want 
            # mismatch values to be returned
            elif identity == 0.0 and comparison[1]+comparison[2] > final_comparison[1]+final_comparison[2]:
                final_comparison = comparison
                final_ref_mRNA = mRNA_ref_id
                final_alt_mRNA = mRNA_alt_id 
    
    # return the highest identity between the mRNA of both locus
    final_identity = round(final_identity * 100, 1)    
    if verbose:
        print(f"\nResult of the comparison of the locus : {final_comparison[1]} matches and {final_comparison[0]} mismatches" )
    return (final_comparison, final_identity, final_ref_mRNA, final_alt_mRNA)
    

## This function compares the loci of the reference and alternative clusters
# returned by the construct_clusters function by assigning each locus to the 
# overlapping locus of the other annotation cluster which gives the highest
# computed identity. Each comparison results are written to a 'results' 
# dictionary detailing locus information.
#
# @see construct_clusters()
#
# @param cluster Instance of class 'Cluster' as returned by construct_clusters
#
# @param create_strings Boolean indicating wether to use the new compare_loci
# fucntion ('False') or the old_compare_loci function ('True')
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @see compare_loci()
#
# @see old_compare_loci()
#
# @returns Returns a list of dictionaries detailing the locus information for 
# each locus comparison 
def annotation_match(cluster, create_strings=False, verbose=False):
    cluster_ref = cluster.get_loci()["ref"]
    cluster_alt = cluster.get_loci()["alt"]
    cluster_name = cluster.name
    if verbose:
        print(f"matching annotations loci with each other for {cluster_name}")
    dyn_prog_matrix = [] # dynamic programmation matrix
    
    # if the loci stored in the cluster are on the reverse strand, reverse 
    # their mRNA's cds lists
    try:
        ref_is_rev = cluster_ref[0].direction == "reverse" 
    except IndexError:
        ref_is_rev = False
    try:
        alt_is_rev = cluster_alt[0].direction == "reverse" 
    except IndexError:
        alt_is_rev = False
    if ref_is_rev == True and alt_is_rev == True:
        cluster_end = cluster.get_end()
        for loc in cluster_ref:
            loc.reverse(cluster_end)
        for loc in cluster_alt:
            loc.reverse(cluster_end)
            
    # initialisation (expand matrix and fill it with zeros)
    for i in range(len(cluster_ref)+1):
        dyn_prog_matrix.append([])
        for j in range(len(cluster_alt)+1):
            dyn_prog_matrix[i].append(0)
            
    # compute all internal values of the matrix by taking the maximum value of
    # the top 'cell', the left 'cell', and the top-left 'cell' summed with the
    # computed identity for its corresponding loci
    for i in range(1, len(cluster_ref)+1):
        for j in range(1, len(cluster_alt)+1):
            
            if create_strings:
                comparison, identity, ref_mRNA, alt_mRNA = old_compare_loci(cluster_ref[i-1], cluster_alt[j-1], verbose)
            else:
                comparison, identity, EI_RF_mismatch_zones, ref_mRNA, alt_mRNA = compare_loci(cluster_ref[i-1], cluster_alt[j-1], verbose)
            
            dyn_prog_matrix[i][j] = max(dyn_prog_matrix[i-1][j],
                                        dyn_prog_matrix[i][j-1],
                                        dyn_prog_matrix[i-1][j-1] + identity)
        
    # retrieve best match alignment through backtracking
    i = len(cluster_ref)
    j = len(cluster_alt)
    results = []
    while i>0 and j>0:
        
        if create_strings:
            comparison, identity, ref_mRNA, alt_mRNA = old_compare_loci(cluster_ref[i-1], cluster_alt[j-1], False)    
        else:
            comparison, identity, mismatch_zones, ref_mRNA, alt_mRNA = compare_loci(cluster_ref[i-1], cluster_alt[j-1], False)
            if cluster_ref[0].direction == "reverse":
                if mismatch_zones not in ["_", "?"]:
                    mismatch_zones = reverse_coord(mismatch_zones, cluster_end)
        
        # if the maximum value is the diagonal value + computed identity, we
        # add both locus informations to the results dictionary and 'progress'
        # to the top-left cell
        if max(dyn_prog_matrix[i-1][j],
            dyn_prog_matrix[i][j-1],
            dyn_prog_matrix[i-1][j-1] + identity) == dyn_prog_matrix[i-1][j-1]+identity:
            if create_strings: 
                if comparison == '_':
                    mismatch_zones = '_'
                else:
                    mismatch_zones = '?'
                
            # we retrieve the number of mRNAs of each locus
            num_mRNAs_ref = len(cluster_ref[i-1].mRNAs.keys())
            num_mRNAs_alt = len(cluster_alt[j-1].mRNAs.keys())
            
            results.append({"reference" : cluster_ref[i-1].name,
                            "reference start" : cluster_ref[i-1].start,
                            "reference end" : cluster_ref[i-1].end,
                            "alternative" : cluster_alt[j-1].name,
                            "alternative start" : cluster_alt[j-1].start,
                            "alternative end" : cluster_alt[j-1].end,
                            "mismatch/match" : comparison,
                            "identity" : identity,
                            "mismatch zones" : mismatch_zones,
                            "cluster name" : cluster_name,
                            "reference mRNA" : ref_mRNA,
                            "alternative mRNA" : alt_mRNA,
                            "reference mRNA number" : num_mRNAs_ref,
                            "alternative mRNA number" : num_mRNAs_alt})
            i -= 1
            j -= 1
        
        # if the maximum value is the top value, we add the reference locus 
        # informations to the results dictionary and 'progress' to the top cell
        elif max(dyn_prog_matrix[i-1][j],
            dyn_prog_matrix[i][j-1],
            dyn_prog_matrix[i-1][j-1] + identity) == dyn_prog_matrix[i-1][j]:
            if create_strings:
                mismatch_zones = '_'
                
            # we retrieve the number of mRNAs of the reference locus
            num_mRNAs_ref = len(cluster_ref[i-1].mRNAs.keys())
            num_mRNAs_alt = '_'
            
            results.append({"reference" : cluster_ref[i-1].name,
                            "reference start" : cluster_ref[i-1].start,
                            "reference end" : cluster_ref[i-1].end,
                            "alternative" : "~",
                            "alternative start" : "_",
                            "alternative end" : "_",
                            "mismatch/match" : [],
                            "identity" : "_",
                            "mismatch zones" : "_",
                            "cluster name" : cluster_name,
                            "reference mRNA" : ref_mRNA,
                            "alternative mRNA" : "_",
                            "reference mRNA number" : num_mRNAs_ref,
                            "alternative mRNA number" : num_mRNAs_alt})
            i -= 1
           
        # if the maximum value is the left value, we add the alternative locus 
        # informations to the results dictionary and 'progress' to the left cell
        else:
            if create_strings:
                mismatch_zones = '_'
                
            # we retrieve the number of mRNAs of the alternative locus
            num_mRNAs_ref = '_'
            num_mRNAs_alt = len(cluster_alt[j-1].mRNAs.keys())
            
            results.append({"reference" : "~",
                            "reference start" : "_",
                            "reference end" : "_",
                            "alternative" : cluster_alt[j-1].name,
                            "alternative start" : cluster_alt[j-1].start,
                            "alternative end" : cluster_alt[j-1].end,
                            "mismatch/match" : [],
                            "identity" : "_",
                            "mismatch zones" : "_",
                            "cluster name" : cluster_name,
                            "reference mRNA" : "_",
                            "alternative mRNA" : alt_mRNA,
                            "reference mRNA number" : num_mRNAs_ref,
                            "alternative mRNA number" : num_mRNAs_alt})
            j -= 1
            
    # if we reached the first line but did not reach the first column,
    # we add the rest of the alternative loci to the results
    while i==0 and j!=0:
    
        # we retrieve the number of mRNAs of the alternative locus
        num_mRNAs_ref = '_'
        num_mRNAs_alt = len(cluster_alt[j-1].mRNAs.keys())

        results.append({"reference" : "~",
                        "reference start" : "_",
                        "reference end" : "_",
                        "alternative" : cluster_alt[j-1].name,
                        "alternative start" : cluster_alt[j-1].start,
                        "alternative end" : cluster_alt[j-1].end,
                        "mismatch/match" : [],
                        "identity" : 0.0,
                        "mismatch zones" : "_",
                        "cluster name" : cluster_name,
                        "reference mRNA" : "_",
                        "alternative mRNA" : "_",
                        "reference mRNA number" : num_mRNAs_ref,
                        "alternative mRNA number" : num_mRNAs_alt})
        j -= 1
        
    # if we reached the first column but did not reach the first line,
    # we add the rest of the reference loci to the results
    while i!=0 and j==0:
    
        # we retrieve the number of mRNAs of the reference locus
        num_mRNAs_ref = len(cluster_ref[i-1].mRNAs.keys())
        num_mRNAs_alt = '_'
        
        results.append({"reference" : cluster_ref[i-1].name,
                        "reference start" : cluster_ref[i-1].start,
                        "reference end" : cluster_ref[i-1].end,
                        "alternative" : "~",
                        "alternative start" : "_",
                        "alternative end" : "_",
                        "mismatch/match" : [],
                        "identity" : 0.0,
                        "mismatch zones" : "_",
                        "cluster name" : cluster_name,
                        "reference mRNA" : "_",
                        "alternative mRNA" : "_",
                        "reference mRNA number" : num_mRNAs_ref,
                        "alternative mRNA number" : num_mRNAs_alt})
        i -= 1    
            
    # sort the results dictionary list so they are in order of locus start
    final_results = sorted(results, key=lambda d: d['reference'])
        
    return final_results
    
    
