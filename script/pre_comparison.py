#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 21:05:26 2024

@author: vetea
"""

import cluster as cl


## Checks if both given annotations dictionaries are sorted by their start
# coordinate in ascending order
#
# @param dict_ref Reference annotation's dictionary, as returned by 
# get_gff_borders
#
# @param dict_alt Alternative annotation's dictionary, as returned by
# get_gff_borders
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @see get_gff_borders()
#
# @returns Returns 'True' if both dictionaries are already sorted, else 'False'
def ref_alt_is_sorted(dict_ref, dict_alt, verbose=False):
    is_sorted = True
    ref_list = list(dict_ref.items())
    alt_list = list(dict_alt.items())
    
    for i in range(len(ref_list)-1):
        if ref_list[i] > ref_list[i+1]:
            is_sorted = False
            
    for j in range(len(alt_list)-1):
        if alt_list[j] > alt_list[j+1]:
            is_sorted = False

    return is_sorted        


## This function creates a list of all the locus coordinates for both 
# annotations and sorts it in ascending order by their lower bound position if
# necessary
#
# @param dict_ref Dictionary containing all loci of the reference annotation, 
# as returned by the 'get_gff_borders' function
#
# @param dict_alt Dictionary containing all loci of the alternative annotation,
# as returned by the 'get_gff_borders' function
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @see get_gff_borders()
#
# @return Returns a list of tuples containing the lower and upper bounds of 
# each locus, its locus ID, and a boolean indicating if the locus was retrieved
# from the reference (True) or the alternative (False)
#
# @remark If the genes detailed in the input GFF files of the program are 
# already sorted in ascendign order, this function's complexity reduces to 
# O(nb_loci), else it is O(nb_loci x log(nb_loci)) (where 'nb_loci' is the
# number of loci in the dictionaries)
def annotation_sort(dict_ref, dict_alt, verbose=False):
    if verbose:
        print("\n\n**************** Constructing the locus order list of the two annotations ****************")
    is_sorted = ref_alt_is_sorted(dict_ref, dict_alt, verbose=False)
    if is_sorted:
        if verbose: print("loci are already sorted")
        locus_order = []
        i = 0
        j = 0
        ref_key_list = list(dict_ref.keys())
        alt_key_list = list(dict_alt.keys())
        
        while i<len(ref_key_list) and j<len(alt_key_list):
            loc_ref = dict_ref[ref_key_list[i]]
            loc_alt = dict_alt[alt_key_list[j]]
            if loc_ref.start < loc_alt.start:
                locus_order.append((loc_ref.start, loc_ref.end, loc_ref.name, True, loc_ref.direction))
                i += 1
            else:
                locus_order.append((loc_alt.start, loc_alt.end, loc_alt.name, False, loc_alt.direction))
                j += 1
                
        if i<len(ref_key_list):
            for k in range(i, len(ref_key_list)):
                loc_ref = dict_ref[ref_key_list[k]]
                locus_order.append((loc_ref.start, loc_ref.end, loc_ref.name, True, loc_ref.direction))
        if j<len(alt_key_list):
            for k in range(j, len(alt_key_list)):
                loc_alt = dict_alt[alt_key_list[k]]
                locus_order.append((loc_alt.start, loc_alt.end, loc_alt.name, False, loc_alt.direction))
    else:
        if verbose: print("loci not sorted, sorting the loci list")
        locus_order = []
        
        # get all reference loci bounds and locus ids as tuples in a list
        # 'True' indicates the tuple is from the reference
        for locus_id, locus in dict_ref.items():
            locus_order.append((locus.start, locus.end, locus.name, True, locus.direction)) 
     		
        # get all alternative loci bounds and locus ids as tuples in a list
        # 'False' indicates the tuple is from the alternative
        for locus_id, locus in dict_alt.items():
            locus_order.append((locus.start, locus.end, locus.name, False, locus.direction)) 
            
        if verbose:
            print("\nSorting the locus order list")
     	
        # sorts the list using the first value of each tuple 
        # (= lower bound of the locus)
        locus_order.sort()
    
    return locus_order


## This function groups every locus in the tuple list 'locus_order' in 
# instances of the 'Clusters' class depending on the overlap of the reference 
# and alternative loci. Each cluster contains a group of mutually-overlapping 
# loci which don't overlap with others.
#
# @see Clusters
#
# @param dict_ref Dictionary of reference loci, as returned by get_gff_borders
#
# @param dict_alt Dictionary of alternative loci, as returned by get_gff_borders
#
# @param locus_order list of tuples containing the information for each
# locus of each annotation in ascending order of their lower border, as 
# returned by annotation_sort
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @see get_gff_borders()
#
# @see annotation_sort()
#
# @remark The clusters names will not necessarily follow each other. 'Cluster6'
# might follow 'Cluster2'
#
# @returns an instance of the Clusters class describing the cluster structure
# of the two annotations loci
def construct_clusters(dict_ref, dict_alt, locus_order, verbose=False):

    # initialisation of the dictionary keeping track of what loci have already
    # been fused
    already_grouped = []
    
    # initialising the return dictionary of instances of 'Cluster' class
    cluster_list = {}
    
    # for each locus in the loci list...
    for i in range(len(locus_order)):   
        # we get the identifier of the current locus, the lower and upper 
        # bounds of the current locus search, and to which annotation it 
        # belongs
        locus_id = locus_order[i][2] 
        locus_borders = [locus_order[i][0], locus_order[i][1]] 
        locus_is_ref = locus_order[i][3]
        locus_strand = locus_order[i][4]
            
        if locus_is_ref:
            is_ref_marker = "_ref"
        else:
            is_ref_marker = "_alt"
            
        # if the current locus has not yet been grouped, create a new cluster,
        # initialize it with empty lists, and add the current locus to it
        if locus_is_ref:
            if locus_id + "_ref" not in already_grouped:
                cluster_name = "cluster " + str(i)
                
                if cluster_name not in cluster_list.keys():
                    cluster = cl.Cluster(cluster_name)
                    cluster_list[cluster_name] = cluster
                cluster.append_to_loci("ref", dict_ref[locus_id])
                already_grouped.append(locus_id + is_ref_marker)
        else:
            if locus_id + "_alt" not in already_grouped:
                cluster_name = "cluster " + str(i)
                
                if cluster_name not in cluster_list.keys():
                    cluster = cl.Cluster(cluster_name)
                    cluster_list[cluster_name] = cluster
                cluster.append_to_loci("alt", dict_alt[locus_id])
                already_grouped.append(locus_id + is_ref_marker)
            
        # while we did not reach the end of the list and a 'stop signal' (j=-1)
        # is not given, for each locus after the current one...
        j = 1
        while j != -1 and j <= len(locus_order)-i-1:
            # we get the identifier, lower bound, and parent annotation of the
            # locus pointed to by the 'j' iterator
            next_locus_id = locus_order[i+j][2]
            next_lower_border = locus_order[i+j][0]
            next_is_ref = locus_order[i+j][3]
            next_strand = locus_order[i+j][4]
            if next_is_ref:
                next_is_ref_marker = "_ref"
            else:
                next_is_ref_marker = "_alt"
                
            # if the locus pointed by 'j' has not already been grouped and 
            # is not on a different strand...
            if next_locus_id + next_is_ref_marker not in already_grouped and locus_strand == next_strand:
                # if the lower bound of the locus pointed by 'j' is inside the
                # current locus' bounds and hasn't yet been included in 
                # a cluster, we add it to the current cluster and to the 
                # 'already_grouped' list
                # upper bound of the current locus search is also extended 
                # to upper bound of the 'j' locus
                if locus_borders[0] <= next_lower_border <= locus_borders[1] or locus_borders[0] >= next_lower_border >= locus_borders[1]:
                    new_upper_border = locus_order[i+j][1]
                    locus_borders[1] = new_upper_border                
                    
                    if next_is_ref:
                        cluster.append_to_loci("ref", dict_ref[next_locus_id])
                        already_grouped.append(next_locus_id + next_is_ref_marker)
                    else:
                        cluster.append_to_loci("alt", dict_alt[next_locus_id])
                        already_grouped.append(next_locus_id + next_is_ref_marker)
                    j += 1
                      
                # if no locus is found in the current locus search's bounds, 
                # a 'stop signal' is given to stop the 'j' loop
                else:
                    j = -1
                    
            # if the locus pointed by 'j' has already been grouped, the upper 
            # bound of the current locus search is extended to upper bound of
            # the 'j' locus 
            elif next_locus_id + next_is_ref_marker in already_grouped:
                new_upper_border = locus_order[i+j][1]
                locus_borders[1] = new_upper_border
                j += 1
                
            # if the locus pointed by 'j' is on a different strand from the 
            # current locus, it is ignored
            elif locus_strand != next_strand:
                j += 1
                         
        # assign the end coordinate value to the cluster
        try:
            ref_end = cluster.get_loci()["ref"][-1].end
        except IndexError:
            ref_end = -1
        try:
            alt_end = cluster.get_loci()["alt"][-1].end
        except IndexError:
            alt_end = -1
        cluster.set_end(max(ref_end, alt_end))
    
    return cluster_list

## Creates a list detailing the codon position of the given annotation 
# ('cds_bounds')at the start of every intersection interval in the given 
# interval list ('area_bounds')
#
# @param cds_bounds List of CDS intervals of an annotation
#
# @param area_bounds List of intersection intervals from two annotations
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @returns Returns a list of numbers describing the codon position of the 
# annotation at the start of the corresponding intersection interval
#
# @remark The intervals must include the upper bounds. For the intersection
# intervals, this can be obtained through the use of the 
# 'get_intervals_with_included_ub' method 
#
# @see get_intervals_with_included_ub()
def get_reading_frame(cds_bounds, area_bounds, verbose=False):
    nb_nt = 0;
    nb_nt_in_cds=0;
    cdsb=0; # cds bounds index
    reading_frames = [];
    for i in range(0, len(area_bounds)-1, 2):
        # move to the CDS containing the current area
        while area_bounds[i] > cds_bounds[cdsb+1]:
            nb_nt += cds_bounds[cdsb+1] - cds_bounds[cdsb] +1;
            cdsb+=2
        # now set the reading frame of the current are
        nb_nt_in_cds = area_bounds[i]-cds_bounds[cdsb]+1;
        reading_frames.append((nb_nt+nb_nt_in_cds-1) % 3 +1);
    return reading_frames;


## This function expects a list of all CDS coordinates (start and end) of a 
# locus. It returns a list indicating the start of the locus as first value 
# and a string describing the codon position of each nucleotide (1,2,3, or 0 
# in the case of a non-CDS nucleotide) of the locus/gene as second value
#
# @param borders The list containing all start-end coordinates of the 
# annotation's locus' CDS
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @see get_gff_borders()
#
# @returns Returns a list describing the start of the locus and the annotation 
# structure of the locus
def create_vectors(borders, verbose=False):
    
    # this variable takes in the strings of gene annotation structure 
    # for each gene
    vector = [0, ""]
    
    if verbose:
        print("\n**************** Converting the coordinates of the locus into a structure string ****************")
    
    # if the locus is on the direct strand
    if borders[1] > borders[0]:
        
        # we get the start of the locus
        vector[0] = borders[0]
        
    # if the locus is on the reverse strand
    else:
        
        # we get the start (the last CDS value since the locus is reversed) 
        # of the locus
        vector[0] = borders[-1]
    
    # this variable takes the codon position of the next CDS nucleotide and 
    # loops between the values 1, 2, and 3
    codon_pos = 1 
    
    # this variable indicates if we are in an exon/CDS or not.
    # it is incremented at each transition between annotations (each element 
    # in the 'borders' list) to represent the exon-intron change along the gene
    in_exon = 1

    # for each coordinate indexed for this locus in 'borders'
    for i in range(len(borders)-1):
    
        # if we are in a CDS, we append the numbers 1, 2, and 3 
        # (with looping) 
        if in_exon % 2 == 1:
                
            # for each nucleotide between this coordinate and the next...
            for j in range( borders[i+1] - borders[i]+1):
                
                # we append the codon position to the structure string 
                vector[1] += str(codon_pos)
                
                # we increment the codon position with looping
                if codon_pos == 3:
                    codon_pos = 1
                else:
                    codon_pos += 1
            
        # if we are not in a CDS, we append 0
        if in_exon % 2 == 0:
            
            # for each nucleotide between this coordinate and the next...
            for j in range( borders[i+1] - borders[i]-1):
                vector[1] += "0"
        
        in_exon += 1    
    
    if verbose:
        print("\nStructure string of the locus :\n" + vector[1] + "\n")
        
    return vector