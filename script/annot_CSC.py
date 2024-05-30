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


## This class represents clusters of overlapping loci computed by the 
# function construct_clusters. It posesses only one attribute, 'clusters',
# which is a dictionary of dictionaries of instances of the class 'Locus'
# indicating the loci attributed to each cluster.
#
# @see construct_clusters()
#
# @see Locus
class Clusters:
    
    ## This method initialises the class with an empty dictionary
    def __init__(self):
        self.clusters = {}
        
    ## This method is used to retrieve the 'clusters' attribute of an instance
    #
    # @returns Returns the 'clusters' attribute of the class instance
    def clusters(self):
        return self.clusters
    
    ## This method is used to retrieve the complete list of mRNAs of the locus
    # of the cluster specified with its locus ID and a boolean indicating
    # if the locus is to be searched for in the reference or alternative list
    #
    # @param loc_id Identifier of the locus from which to retrieve mRNAs
    #
    # @param ref Boolean indicating from which of the reference (True) or 
    # alternative (False) locus to retrieve the mRNAs
    #
    # @returns Returns the list of all mRNAs of the specified locus
    #
    # @remark This method is intended to be used in unit tests to verify
    # the clusters contains the intended loci
    def get_mRNAs(self, loc_id, ref):
        list_mRNAs = []
        for loc in self.clusters[loc_id][ref]:
            list_mRNAs.append(loc.mRNAs)
        return list_mRNAs
    

## This class represents an annotation's locus identified from reading the
# associated GFF file with the get_gff_borders function. It posesses multiple 
# attributes describing the locus
#
# @see get_gff_borders()
class Locus:
    
    ## This method initialises the class with default values for each attribute
    #
    # @param name ID of the locus
    #
    # @param mRNAs Dictionary listing all mRNAs of the locus. Each key 
    # corresponds to the ID of the mRNA and eahc value to the list of all CDS
    # coordinates as retrieved by the get_gff_borders function
    #
    # @param start Start coordinate of the locus
    #
    # @param end End coordinate of the locus
    #
    # @param dicrection String indicating if the locus is on the 'direct' or 
    # 'reverse' strand
    #
    # @see get_gff_borders()
    def __init__(self, name="", mRNAs=None, start=-1, end=-1, direction=""):
        if not mRNAs:
            mRNAs = {}
        
        self.name = name
        self.mRNAs = mRNAs.copy()
        self.start = start
        self.end = end
        self.direction = direction
        
    ## This method is used to retrieve the 'mRNAs' attribute of an instance
    #
    # @returns Returns the 'mRNAs' attribute of the class instance
    def mRNAs(self):
        return self.mRNAs
    
    ## This method is used to verify if all mRNAs of a dictionary are present 
    # in the class instance's 'mRNAs' attribute
    #
    # @param **mRNAs Dictionary with mRNA names as keys and CDS coordinates
    # list as values
    #
    # @returns Returns a boolean indicating if all given mRNAs were found in
    # the class instance's 'mRNAs' attribute
    #
    # @remark This method is intended to be used in unit tests to verify
    # the locus contains the intended mRNAs
    def contain_mrnas(self, **mrnas):
        for mrna_name, positions_list in mrnas.items():
            if mrna_name not in self.mRNAs:
                return False
            return self.mRNAs[mrna_name] == positions_list

    ## This function was added for convenience as a way to easily retrieve the
    # attribute values of the class instance as a formatted string
    #
    # @returns Returns a string describing the values of each attribute of the
    # class instance
    def show_init(self):
        return f"Locus(name='{self.name}', mRNAs={self.mRNAs}, start={self.start}, end={self.end}, direction='{self.direction}'"


## This function retrieves and returns the id of the structure described from 
# a line read from a GFF file.
#
# @param line The line read from the file (string)
#
# @param debug If True, triggers display of many messages intended for 
# debugging the program. Default is 'False'
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @remark This function expects the file to be in GFF format
#
# @returns Returns the id of the structure described by the line
def get_structure_id(line, debug=False, verbose=False):
    
    # try to retrieve the last column, else return an error
    try: 
        last_col = line.split("\t")[8]
    except IndexError:
        print("\nFile is not in GFF format (no ninth column)")
        sys.exit(1)
    
    # retrieve the id with the rest of the column text (commentaries)
    try:
        id_and_rest = last_col.split("=")[1]
    except IndexError:
        print("\nNo ID field found in last column (ninth column)")
        sys.exit(1)

    structure_id = ""
    i = 0
    
    # for each character from the first, we add it to the locus_id if it is 
    # not in a list of special characters. When the first special character is 
    # encountered, stop the loop and return the locus_id 
    while id_and_rest[i] not in [",", "?", ";", ":", "/", "!", "*", "$", "%", "+", "@", "#", "~", "&", "\n", "\t"] :
        structure_id += id_and_rest[i]
        i += 1
    if debug:
        print(f"Structure ID = {structure_id}")
        
    return structure_id


## This function retrieves and returns the parent id of the structure 
# described from a line read from a GFF file.
#
# @param line The line read from the file (string)
#
# @param debug If True, triggers display of many messages intended for 
# debugging the program. Default is 'False'
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @remark This function expects the file to be in GFF format, and the 
# parent id field to come after the id field
#
# @returns Returns the id of the parent of the structure described by the line
def get_parent_id(line, debug=False, verbose=False):
    
    # try to retrieve the last column, else return an error
    try: 
        last_col = line.split("\t")[8]
    except IndexError:
        print("\nFile is not in GFF format (no ninth column)")
        sys.exit(1)
    
    # retrieve the id with the rest of the column text (commentaries)
    try:
        id_and_rest = last_col.split("=")[2]
    except IndexError:
        print("\nNo parent ID field found in last column (ninth column)")
        sys.exit(1)
    
    structure_id = ""
    i = 0
    
    # for each character from the first, we add it to the locus_id if it is 
    # not in a list of special characters. When the first special character is
    # encountered, stop the loop and return the locus_id 
    while id_and_rest[i] not in [",", "?", ";", ":", "/", "!", "*", "$", "%", "+", "@", "#", "~", "&", "\n", "\t"] :
        if debug:
            print(f"Reading character {id_and_rest[i]}")
        structure_id += id_and_rest[i]
        i += 1
    if debug:
        print(f"Structure parent ID = {structure_id}")
        
    return structure_id
    


## This function expects a string corresponding to the file path of the GFF 
# file to read, and returns a dictionary of instances of the class 'Locus', 
# detailing all the relevant information for the gene and its mRNAs
#
# @param path Path of the file to read
#
# @param debug If True, triggers display of many messages intended for 
# debugging the program. Default is 'False'
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @remark If the parent ID of a CDS does not match the ID of the previous 
# mRNA (indicating an incorrect file structure), an entry is added to a 'log' 
# file but the function is not interrupted
#
# @return Returns a dictionary of instances of the class 'Locus', containing 
# the information of the CDS borders of each mRNA of the gene, the start and
# end coordinates, the DNA strand on which the gene is predicted, and the
# locus ID
#
# @see Locus
#
# @remark if a locus is on the reverse strand, its CDS borders list is reversed
# to be in ascending order
#
def get_gff_borders(path, debug=False, verbose=False):
    
    try:    
        file = open(path, "r") # the file to read
    except FileNotFoundError:
        print(f"\nget_gff_borders() function error : file '{path}' does not exist or cannot be read from\n")
        sys.exit(2)
    
    # try to open the structure errors log file
    try:
        log = open("./results/log.txt", "w") 
    except FileNotFoundError:
        os.mkdir("./results/") # create 'results' subdirectory
        log = open("./results/log.txt", "w")
        
    loci = {} # return dictionary
    locus_id = "" # the current gene being analyzed
    mRNA_id = "" # the current mRNA being analyzed
    line_index = 1 # number of the file line currently read
    locus = Locus() # Initialisation of the Locus class instance to construct
    
    for line in file:
        
        # if we encounter a new gene which is not from a contig, 
        # we get its information and create an instance in 'loci'
        if str(line.split("\t")[2]) == "gene" and line.split("\t")[0] != "contig": 
            
            # if there was a previous gene, but its borders list is empty, 
            # return an error
            if locus_id != "" and locus.mRNAs[mRNA_id] == []:
                print(f"\nLine {line_index} = get_gff_borders() function error : no coding sequence (CDS) could be found for the previous mRNA '{mRNA_id}'\n")
                sys.exit(1)
            
            # if there was a previous gene, add it to return dictionary and
            # delete it. Then we create a new Locus instance
            # (if the gene is on the reverse strand, reverse its border list
            # before adding it to the dictionary)
            if locus.mRNAs != {}:
                if locus.direction == "reverse":
                    locus.mRNAs[mRNA_id].reverse()
                loci[locus_id] = locus
            del locus
            locus = Locus()
            
            # we retrieve the id of the locus
            locus_id = get_structure_id(line, debug, verbose)
            locus.name = locus_id
            
            # we retrieve the start and end coordinates of the locus
            start = int(line.split("\t")[3])
            locus.start = start
            end = int(line.split("\t")[4])
            locus.end = end
            
            # we deduce from the coordinates the direction of the locus
            if locus.end < locus.start:
                locus.direction = "reverse"
            else:
                locus.direction = "direct"
            
            if verbose :
                print("\n**************** Reading the locus " + locus_id + " ****************")

        # if we encounter a new mRNA, we get its ID and create a key 
        # in the instance's mRNAs attribute with an empty list
        if str(line.split("\t")[2]) == "mRNA" and line.split("\t")[0] != "contig": 
            
            mRNA_id = get_structure_id(line, debug, verbose)
            locus.mRNAs[mRNA_id] = []
            if verbose :
                print("\n**************** Reading mRNA " + locus_id + " ****************")

        # if we encounter a CDS, we add its start and end positions to the 
        # corresponding mRNA key in the locus' mRNAs attribute
        if str(line.split("\t")[2]) == "CDS" and line.split("\t")[0] != "contig":
            parent_id = get_parent_id(line, debug, verbose)
            
            if parent_id != mRNA_id:
                print("\nIncorrect file structure (Parent of CDS is not previous mRNA). See 'log.txt' for more information")
                log.write("Line " + str(line_index) + " : CDS parent ID (" + parent_id + ") does not match last mRNA ID (" + locus_id +")\n")
                
            if locus_id == '':
                print(f"\nLine {line_index} = get_gff_borders() function error : CDS has been found before any mRNA")
                sys.exit(1)
            
            locus.mRNAs[mRNA_id].append(int(line.split("\t")[3]))
            locus.mRNAs[mRNA_id].append(int(line.split("\t")[4]))
            if verbose:
                print("Adding borders to " + locus_id + " : " + line.split("\t")[3] + ", " + line.split("\t")[4])
                
        line_index += 1 # increment the line indicator (for debug purposes)
                
    if locus_id == "": # if the locus_id still has default value, return error
        print(f"\nget_gff_borders() function error : no locus could be found in file '{path}'\n")
        sys.exit(1)
        
    file.close()
    log.close()
    if locus.direction == "reverse":
        locus.mRNAs[mRNA_id].reverse()
    loci[locus_id] = locus # add the last locus to the return dictionary
    
    # we return the entire dictionary with all loci
    return loci


## This function creates a list of all the locus coordinates for both 
# annotations and sorts it in ascending order by their lower bound position. 
#
# @param dict_ref Dictionary containing all loci of the reference annotation, 
# as returned by the 'get_gff_borders' function
#
# @param dict_alt Dictionary containing all loci of the alternative annotation,
# as returned by the 'get_gff_borders' function
#
# @param debug If True, triggers display of many messages intended for 
# debugging the program. Default is 'False'
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @see get_gff_borders()
#
# @return Returns a list of tuples containing the lower and upper bounds of 
# each locus, its locus ID, and a boolean indicating if the locus was retrieved
# from the reference (True) or the alternative (False)
def annotation_sort(dict_ref, dict_alt, debug=False, verbose=False):
    if verbose:
        print("\n**************** Constructing the locus order list of the two annotations ****************")
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
    if debug:
        print(f"\nlocus order list = {locus_order}")
     
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
# @param debug If True, triggers display of many messages intended for 
# debugging the program. Default is 'False'
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
# @remark Loci on different strands ('direct' or 'reverse') will be put in
# different clusters
#
# @returns an instance of the Clusters class describing the cluster structure
# of the two annotations loci
def construct_clusters(dict_ref, dict_alt, locus_order, debug=False, verbose=False):

    # initialisation of the dictionary keeping track of what loci have already
    # been fused
    already_grouped = []
    
    # initialising a new instance of 'Clusters' class
    clusters = Clusters()
    
    # for each locus in the loci list...
    for i in range(len(locus_order)):   
        # we get the identifier of the current locus, the lower and upper 
        # bounds of the current locus search, and to which annotation it 
        # belongs
        locus_id = locus_order[i][2] 
        locus_borders = [locus_order[i][0], locus_order[i][1]] 
        locus_is_ref = locus_order[i][3]
        locus_strand = locus_order[i][4]
        if debug:
            print(f"\ni = {i}")
            print(f"locus id = {locus_id}")
            print(f"locus borders = {locus_borders}")
            print(f"locus strand = {locus_strand}")
            
        if locus_is_ref:
            is_ref_marker = "_ref"
        else:
            is_ref_marker = "_alt"
            
        # if the current locus has not yet been grouped, create a new cluster,
        # initialize it with empty lists, and add the current locus to it
        if locus_is_ref:
            if locus_id + "_ref" not in already_grouped:
                cluster_name = "cluster " + str(i)
                
                if cluster_name not in clusters.clusters:
                    clusters.clusters[cluster_name] = {"ref" : [],
                                                   "alt" : []}
                clusters.clusters[cluster_name]["ref"].append(dict_ref[locus_id])
                already_grouped.append(locus_id + is_ref_marker)
        else:
            if locus_id + "_alt" not in already_grouped:
                cluster_name = "cluster " + str(i)
                
                if cluster_name not in clusters.clusters:
                    clusters.clusters[cluster_name] = {"ref" : [],
                                                   "alt" : []}
                clusters.clusters[cluster_name]["alt"].append(dict_alt[locus_id])
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
            if debug:
                print(f"\nj = {j}")
                print(f"next lower border = {next_lower_border}")
                print(f"already grouped : {already_grouped}")
                print(f"next locus strand = {next_strand}")
                
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
                    if debug:
                        print("next locus found in current borders")
                    new_upper_border = locus_order[i+j][1]
                    locus_borders[1] = new_upper_border                
                    if debug:
                        print(f"new locus borders = {locus_borders}")
                    
                    if next_is_ref:
                        clusters.clusters[cluster_name]["ref"].append(dict_ref[next_locus_id])
                        already_grouped.append(next_locus_id + next_is_ref_marker)
                    else:
                        clusters.clusters[cluster_name]["alt"].append(dict_alt[next_locus_id])
                        already_grouped.append(next_locus_id + next_is_ref_marker)
                    j += 1
                      
                # if no locus is found in the current locus search's bounds, 
                # a 'stop signal' is given to stop the 'j' loop
                else: 
                    if debug:
                        print("next locus not found in current borders, stopping search")
                    j = -1
                    
            # if the locus pointed by 'j' has already been grouped, the upper 
            # bound of the current locus search is extended to upper bound of
            # the 'j' locus 
            elif next_locus_id + next_is_ref_marker in already_grouped:
                if debug:
                    print(f"{next_locus_id} already grouped, skipping and extending search")
                new_upper_border = locus_order[i+j][1]
                locus_borders[1] = new_upper_border               
                if debug:
                    print(f"new locus borders = {locus_borders}")               
                j += 1
                
            # if the locus pointed by 'j' is on a different strand from the 
            # current locus, it is ignored
            elif locus_strand != next_strand:
                if debug:
                    print(f"{next_locus_id} on different strand, skipping")        
                j += 1
                         
        if debug:
            print(f"cluster name = {cluster_name}")
            print(clusters.clusters[cluster_name]["ref"])
            print(clusters.clusters[cluster_name]["alt"])
    
    return clusters


## This function retrieves all the CDS coordinates from the given lists of 
# coordinates of the reference and alternative annotations
# and includes them in a unique list of coordinates. The coordinates are 
# sorted in ascending order.
#
# @param ref List of CDS coordinates of the reference annotation
#
# @param alt List of CDS coordinates of the alternative annotation
#
# @param debug If True, triggers display of many messages intended for 
# debugging the program. Default is 'False'
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @return Returns a list of coordinates compiling all coordinates from both 
# initial lists in ascending order
def get_area_bounds(ref, alt, debug=False, verbose=False):
    
    if verbose:
        print("\n**************** Creating comparison areas for the two annotations ****************")
    bounds=[] # the return list
    i=0
    j=0    
    
    # while we did not reach the end of the coordinates lists...
    while i <= len(ref)-1 and j <= len(alt)-1:
        if debug:
            print(f"i = {i}, j = {j}")
            print(f"minimum of ref[i] ({ref[i]}) and alt[j] ({alt[j]}) = {min(ref[i], alt[j])}")
        
        # we add the next coordinate to the result list and increment the 
        # corresponding iterator
        bounds.append(min(ref[i],alt[j]))
        
        if(ref[i] == bounds[-1]):
            i+=1
        if(alt[j] == bounds[-1]):
            j+=1
            
    # after we get to the end of one annotation, we append the rest of the other one to the result list
    
    if i > len(ref)-1:
        for k in range(j, len(alt)):
            bounds.append(alt[k])
        
    if j > len(alt)-1:
        for k in range(i, len(ref)):
            bounds.append(ref[k])
            
    if debug:
        print(f"Final area_bounds = {bounds}")
            
    return bounds

## This function indicates for each couple of bounds in the given list of area
# bounds (@param area_bounds), if they delimit an area which includes a CDS 
# from the given CDS coordinates list (@param cds_bounds). It is used during 
# the comparison of areas in compare_loci() to know if the reference or 
# alternative have a CDS in the area
#
# @see compare_loci()
#
# @param cds_bounds List of CDS coordinates for an annotation
#
# @param area_bounds List of area bounds
#
# @param debug If True, triggers display of many messages intended for 
# debugging the program. Default is 'False'
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @return Returns a list indicating for each couple of bounds if they include 
# a CDS ('True') or not ('False')
def is_in_cds(cds_bounds, area_bounds, debug=False, verbose=False):
    cds_id=0
    i=0
    in_cds=[] # return list
    
    # while we did not yet reach the end of the bounds list or the CDS coordinates list...
    while (i < len(area_bounds)-1 and cds_id < len(cds_bounds)-1):
        
        current_cds = [cds_bounds[cds_id], cds_bounds[cds_id+1]]
        lb = area_bounds[i]
        ub = area_bounds[i+1]
        if debug:
            print(f"\ni = {i}")
            print(f"cds_id = {cds_id}")
            print(f"current_cds = {current_cds}")
            print(f"Lower bound (lb) = {lb}")
            print(f"Upper bound (ub) = {ub}")

        # if the current CDS is included in the area bounds, 'True' is used to index it in the return list, else 'False'
        if(current_cds[0] <= lb and ub <= current_cds[1]):
            in_cds.append(True)
        else:
            in_cds.append(False)
        if debug:
            print(f"in_cds = {in_cds}")
            
        # if the current area reaches the end of the current CDS, we skip to the next CDS
        if(ub == current_cds[1]):
            cds_id += 2
            
        i += 1
        
    # if we reached the end of the CDS coordinates list, we continue until all area bounds are exhausted by appending 'False'
    while(i < len(area_bounds)-1):
        if debug:
            print(f"i (after end of CDS list) = {i}")
        
        in_cds.append(False)
        i += 1
        
    if debug:
        print(f"\nFinal in_cds = {in_cds}\n")
        
    return in_cds


## This function compares two annotations' loci returned by the function 
# get_gff_borders and creates for each pair of reference-alternative mRNAs 
# a comparison list detailing the identities and differences between the two
# annotations's codon position structure. It returns a tuple containing the 
# comparison list giving the highest identity, the computed identity level, 
# and the list of mismatch areas indentified by the comparison.
#
# @see get_gff_borders()
#
# @param ref_locus The reference annotation's locus (Locus class instance)
#
# @param alt_locus The alternative annotation's locus (Locus class instance)
#
# @param debug If True, triggers display of many messages intended for 
# debugging the program. Default is 'False'
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @returns Returns a tuple containing a list indicating the number of codon 
# position mismatches (first position) and matches (second position) between 
# the two border lists giving the maximum indentity between all mRNAs, the 
# identity level computed from this list, and the comaprison areas
# producing mismatches
#
# @remark This function doesn't expect any annotation to be a 'reference'
def compare_loci(ref_locus, alt_locus, debug=False, verbose=False):
    
    if verbose:
        print(f"\n**************** comparing loci {ref_locus.name} of reference and {alt_locus.name} of alternative ****************")
        
    # initialize final return values
    final_comparison = [0,0]
    final_identity = 0.0
    final_mismatch_zones = []
    
    # if the loci are on different strands, return 'null' values
    if ref_locus.direction != alt_locus.direction:
        return ('_', 0.0, '_')
    
    # for each mRNA in the reference locus...
    for mRNA_ref_id, mRNA_ref in ref_locus.mRNAs.items():
        # for each mRNA in the alternative locus...
        for mRNA_alt_id, mRNA_alt in alt_locus.mRNAs.items():
            if debug:
                print(f"mRNA_ref name = {mRNA_ref_id}")
                print(f"mRNA_alt_ name = {mRNA_alt_id}")
                print(f"mRNA_ref = {mRNA_ref}")
                print(f"mRNA_alt = {mRNA_alt}")
                
            if verbose:
                print(f"\n**************** comparing mRNA {mRNA_ref_id} of reference locus {ref_locus.name} and mRNA {mRNA_alt_id} of alternative locus {alt_locus.name} ****************")
            
            # we retrieve the bounds of all areas delimited by all the CDS 
            # coordinates of both border lists
            area_bounds = get_area_bounds(mRNA_ref, mRNA_alt, debug, verbose)
            
            if debug:
                print(f"area bounds = {area_bounds}")
            
            # we retrieve the list indicating the areas which include or not 
            # a CDS for both annotations
            if verbose:
                print("\n**************** Evaluating presence of CDS in the comparison areas ****************")
            ref_in_CDS = is_in_cds(mRNA_ref, area_bounds, debug, verbose)
            alt_in_CDS = is_in_cds(mRNA_alt, area_bounds, debug, verbose)
            if debug:
                print("Ref_in_CDS = " + str(ref_in_CDS))
                print("Alt_in_CDS = " + str(alt_in_CDS))
            
            # initialize the return values for this mRNA
            # (first value of 'comparison' is mismatches, second value 
            # is identities)
            codon_position_ref=0
            codon_position_alt=0
            comparison = [0,0] 
            mismatch_zones = []
            bound_id = 0
            prev_bounds = area_bounds[0]
            
            # for each comparison area delimited by the bounds in the 
            # list 'area_bounds'...
            for bound in area_bounds[1:] :
                if debug:
                    print("\n")
                    print("Bound = " + str(bound))
                    print("Prev_bound = " + str(prev_bounds))
                    print("Bound_id = " + str(bound_id))
                    print("Codon_position_ref = " + str(codon_position_ref))
                    print("Codon_position_alt = " + str(codon_position_alt))
                    print("Ref_in_CDS[bound_id] = " + str(ref_in_CDS[bound_id]))
                    print("Alt_in_CDS[bound_id] = " + str(alt_in_CDS[bound_id]))
                    print("\n")
                    
                # if both annotations are in a CDS and have the same codon 
                # position, then all codon positions for the rest of the area 
                # will be identical, so we add the length of the area to the 
                # second value of the comparison list and update both 
                # codon positions
                if(ref_in_CDS[bound_id] and alt_in_CDS[bound_id] and codon_position_alt == codon_position_ref):
                    if debug:
                        print(f"Identical codon positions for the area, adding {bound-prev_bounds} to match values")
                    
                    comparison[1] += bound-prev_bounds
                    codon_position_alt = (codon_position_alt + (bound-prev_bounds))%3 
                    codon_position_ref = (codon_position_ref + (bound-prev_bounds))%3
                    
                # if both annotations are in a CDS but don't have the same 
                # codon position, then all codon positions for the rest of the
                # area will be different, so we add the length of the are to 
                # the first value of the comparison list and update both codon
                # positions
                elif(ref_in_CDS[bound_id] and alt_in_CDS[bound_id] and codon_position_alt != codon_position_ref):
                    if debug:
                        print(f"Different codon positions for the area, adding {bound-prev_bounds} to mismatch values")
                        
                    comparison[0] += bound-prev_bounds
                    codon_position_alt = (codon_position_alt + (bound-prev_bounds))%3 
                    codon_position_ref = (codon_position_ref + (bound-prev_bounds))%3 
                    mismatch_zones.append(f"{prev_bounds}-{bound}")
                
                # if only one annotation has a CDS in the comparison area, we 
                # add to the first value of the comparison list and update only
                # one codon position
                
                elif(ref_in_CDS[bound_id]):
                    if debug:
                        print(f"Alternative is not in CDS for the area, adding {bound-prev_bounds} to mismatch values")
                        
                    comparison[0] += bound-prev_bounds
                    codon_position_ref = (codon_position_ref + (bound-prev_bounds))%3 
                    mismatch_zones.append(f"{prev_bounds}-{bound}")
                    
                
                elif(alt_in_CDS[bound_id]):
                    if debug:
                        print(f"Reference is not in CDS for the area, adding {bound-prev_bounds} to mismatch values")
                        
                    comparison[0] += bound-prev_bounds
                    codon_position_alt = (codon_position_alt + (bound-prev_bounds))%3 
                    mismatch_zones.append(f"{prev_bounds}-{bound}")
                        
                # the case of both annotations being outside of a CDS is not 
                # used in the computation of global loci identity, and is 
                # ignored
                else:
                    if debug:
                        print("Reference and alternative are not in CDS for the area, ignoring area")
                    
                prev_bounds = bound
                bound_id += 1
                
            # when the loci are on the reverse strand, matches and mismatches
            # are counted in the negatives, so we convert them using the 
            # absolute value
            comparison[0] = abs(comparison[0])
            comparison[1] = abs(comparison[1])    
            
            if verbose:
                print(f"\nResult of the comparison of the locus : {comparison[1]} matches and {comparison[0]} mismatches\n")
                
            identity = comparison[1] / (comparison[0] + comparison[1])
        
            # for each mRNA, we test wether the computed identity is higher 
            # than for the preceding mRNAs, to retrieve the highest identity
            if identity > final_identity:
                final_comparison = comparison
                final_identity = identity
                final_mismatch_zones = mismatch_zones
                if debug:
                    print(f"comparison = {comparison}\nfinal_comparison = {final_comparison}\nidentity = {identity}\nfinal identity = {final_identity}\nmismatch zones = {mismatch_zones}\nfinal mismatch zones = {final_mismatch_zones}")
            # if all mRNAs comparisons return 0% identity, we still want 
            # mismatch values to be returned
            elif identity == 0.0 and comparison[0] > final_comparison[0]:
                final_comparison = comparison
                final_mismatch_zones = mismatch_zones
    
    # return the highest mRNA identity between the locus of each annotation
    # (as a percentage)
    final_identity = round(final_identity * 100, 1)
    if debug:
        print(f"final mismatch zones = {final_mismatch_zones}")
    if verbose:
        print(f"\nFinal result of the comparison of the locus : {final_comparison[1]} matches and {final_comparison[0]} mismatches (identity = {final_identity})" )
    return (final_comparison, final_identity, final_mismatch_zones)


## This function expects a list of all CDS coordinates (start and end) of a 
# locus. It returns a list indicating the start of the locus as first value 
# and a string describing the codon position of each nucleotide (1,2,3, or 0 
# in the case of a non-CDS nucleotide) of the locus/gene as second value
#
# @param borders The list containing all start-end coordinates of the 
# annotation's locus' CDS
#
# @param debug If True, triggers display of many messages intended for 
# debugging the program. Default is 'False'
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @see get_gff_borders()
#
# @returns Returns a list describing the start of the locus and the annotation 
# structure of the locus
def create_vectors(borders, debug=False, verbose=False):
    
    # this variable takes in the strings of gene annotation structure 
    # for each gene
    vector = [0, ""]
    
    if verbose:
        print("\n**************** Converting the coordinates of the locus into a structure string ****************")
    
    # if the locus is on the direct strand
    if borders[1] > borders[0]:
        if debug:
            print(f"\nLocus is on direct strand, retrieving start of locus : {borders[0]}")
        
        # we get the start of the locus
        vector[0] = borders[0]
        
    # if the locus is on the reverse strand
    else:
        if debug:
            print(f"\nLocus is on reverse strand, retrieving start of locus from end of coordinates list : {borders[-1]}")
        
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
    
    if borders[1] > borders[0]:
    
        # for each coordinate indexed for this locus in 'borders'
        for i in range(len(borders)-1):
        
            # if we are in a CDS, we append the numbers 1, 2, and 3 
            # (with looping) 
            if in_exon % 2 == 1:
                if debug:
                    print(f"i = {i}")
                    print("in_exon = True (adding codon positions to structure string)")
                    print(f"Codon position = {codon_pos}")
                    
                # for each nucleotide between this coordinate and the next...
                for j in range( borders[i+1] - borders[i] ):
                    
                    # we append the codon position to the structure string 
                    vector[1] += str(codon_pos)
                    
                    # we increment the codon position with looping
                    if codon_pos == 3:
                        codon_pos = 1
                    else:
                        codon_pos += 1
                if debug:
                    print(f"New codon position = {codon_pos}")
                    print(f"New structure string : {vector[1]}")
                
            # if we are not in a CDS, we append 0
            if in_exon % 2 == 0:
                if debug:
                    print(f"i = {i}")
                    print("in_exon = False (adding 0 to structure string)")
                
                # for each nucleotide between this coordinate and the next...
                for j in range( borders[i+1] - borders[i] ):
                    vector[1] += "0"
                if debug:
                    print(f"New structure string : {vector[1]}")
            
            in_exon += 1    
        
    elif borders[1] < borders[0]:

        # for each coordinate indexed for this locus in 'borders' in 
        # reverse order
        for i in range(len(borders)-1, 0, -1):
        
            # if we are in a CDS, we append the numbers 1, 2, and 3 
            # (with looping) 
            if in_exon % 2 == 1:
                if debug:
                    print(f"i = {i}")
                    print("in_exon = True (adding codon positions to structure string)")
                    print(f"Codon position = {codon_pos}")
                
                # for each nucleotide between this coordinate and the next...
                for j in range( borders[i-1] - borders[i] ):
                    
                    # we append the codon position to the structure string 
                    vector[1] += str(codon_pos)
                    
                    # we increment the codon position with looping
                    if codon_pos == 3:
                        codon_pos = 1
                    else:
                        codon_pos += 1
                if debug:
                    print(f"New codon position = {codon_pos}")
                    print(f"New structure string : {vector[1]}")
                
            # if we are not in a CDS, we append 0
            if in_exon % 2 == 0:
                if debug:
                    print(f"i = {i}")
                    print("in_exon = False (adding 0 to structure string)")
                
                # for each nucleotide between this coordinate and the next...
                for j in range( borders[i-1] - borders[i] ):
                    vector[1] += "0"
                if debug:
                    print(f"New structure string : {vector[1]}")
            
            in_exon += 1   
    
    if verbose:
        print("\nStructure string of the locus :\n" + vector[1] + "\n")
        
    return vector


## This function expects two loci corresponding to two annotations of the same 
# genome, creates and compares their structure strings (create_vectors function) 
#
# @see create_vectors()
#
# @param borders_vector_ref Instance of class 'Locus' of the reference locus
#
# @param borders_vector_alt Instance of class 'Locus' of the alternative locus
#
# @param debug If True, triggers display of many messages intended for 
# debugging the program. Default is 'False'
#
# @param verbose If True, triggers display of more information messages. 
# Default is 'False'
#
# @returns Returns a tuple containing the list of mismatches (first value) 
# and matches (second value) and the computed string identity
#
# @remark This function doesn't expect any annotation to be a 'reference'
def old_compare_loci(ref_locus, alt_locus, debug=False, verbose=False): 
    
    if verbose:
        print(f"\n**************** comparing loci {ref_locus.name} of reference and {alt_locus.name} of alternative ****************")
    
    # initialisation of the final return values
    final_comparison = [0,0]
    final_identity = 0.0
    
    # if the loci are on different strands, return 'null' values
    if ref_locus.direction != alt_locus.direction:
        return ('_', 0.0)
    
    # for each reference locus mRNA...
    for mRNA_ref_id, mRNA_ref in ref_locus.mRNAs.items():
        # for each alternative locus mRNA...
        for mRNA_alt_id, mRNA_alt in alt_locus.mRNAs.items():
                
            if verbose:
                print(f"\n**************** comparing mRNA {mRNA_ref_id} of reference locus {ref_locus.name} and mRNA {mRNA_alt_id} of alternative locus {alt_locus.name} ****************")
            
            # we create the structure string for each mRNA
            start_ref, vector_ref = create_vectors(mRNA_ref, debug, verbose)
            start_alt, vector_alt = create_vectors(mRNA_alt, debug, verbose)
            comparison = [0,0]
                
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
                
            # we get the difference between the start position of each locus
            diff = start_max - start_min
            
            if debug:
                print(f"vector_ref = {vector_ref}")
                print(f"vector_alt = {vector_alt}")
                print(f"vector_min = {vector_min}")
                print(f"vector_max = {vector_max}")
                print(f"position difference = {diff}")
                print(f"range = {max(len(vector_ref), len(vector_alt)) + diff}")
            
            # for each position in the string of maximum length + position diff
            for i in range(max(len(vector_ref), len(vector_alt)) + diff):
                if debug:
                    print(f"i = {i}")
                    
                # identitify the presence of a CDS at each position in
                # the minimum and maximum strings (accounting for position
                # difference)
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
                    if debug:
                        print("both ref and alt are in cds with same codon position value, adding to 'match'")
                    comparison[1] += 1
                    
                # if both are in a cds but don't have the same codon
                # position, identify as a 'mismatch'
                elif min_in_cds and max_in_cds and vector_min[i] != vector_max[i-diff]:
                    if debug:
                        print("both ref and alt are in cds, with different codon position value, adding to 'mismatch'")
                    comparison[0] += 1
                    
                # if only one is in a cds, identify as a 'mismatch'
                elif min_in_cds or max_in_cds:
                    if debug:
                        print("ref or alt in cds, adding to 'mismatch'")
                    comparison[0] += 1
                    
                # else both are not in cds, this case is ignored
                else:
                    if debug:
                        print("ref and alt not in cds, ignoring")
                
            identity = comparison[1] / (comparison[0] + comparison[1])
            
            # for each mRNA, we test wether the computed identity is higher 
            # than for the preceding mRNAs, to retrieve the highest identity
            if identity > final_identity:
                final_comparison = comparison
                final_identity = identity
                if debug:
                    print(f"comparison = {comparison}\nfinal_comparison = {final_comparison}\nidentity = {identity}\nfinal identity = {final_identity}")
            # if all mRNAs comparisons return 0% identity, we still want 
            # mismatch values to be returned
            elif comparison[0] > final_comparison[0]:
                final_comparison = comparison
    
    # return the highest identity between the mRNA of both locus
    # (as a percentage)
    final_identity = round(final_identity * 100, 1)    
    if verbose:
        print(f"\nResult of the comparison of the locus : {final_comparison[1]} matches and {final_comparison[0]} mismatches" )
    return (final_comparison, final_identity)
    

## This function compares the loci of the reference and alternative clusters
# returned by the construct_clusters function by assigning each locus to the 
# overlapping locus of the other annotation cluster which gives the highest
# computed identity. Each comparison results are displayed in the terminal and 
# written to a 'results' dictionary detailing locus information.
#
# @see construct_clusters()
#
# @param cluster_ref List of Locus class instances describing the loci
# of the reference annotation cluster
#
# @param cluster_alt List of Locus class instances describing the loci
# of the alternative annotation cluster
#
# @param cluster_name name of the cluster containing the loci
#
# @param create_strings Boolean indicating wether to use the new compare_loci
# fucntion ('False') or the old_compare_loci function ('True')
#
# @param debug If True, triggers display of many messages intended for 
# debugging the program. Default is 'False'
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
def annotation_match(cluster_ref, cluster_alt, cluster_name, create_strings=False, debug=False, verbose=False):
    if verbose:
        print("\n**************** matching annotations loci with each other ****************")
    dyn_prog_matrix = [] # dynamic programmation matrix
    
    # initialisation (expand matrix and fill it with zeros)
    for i in range(len(cluster_ref)+1):
        dyn_prog_matrix.append([])
        
        for j in range(len(cluster_alt)+1):
            dyn_prog_matrix[i].append(0)
            
    if debug:        
        print("Initialized dynamic programmation matrix :")
        for line in dyn_prog_matrix:
            print(str(line))
            
    # compute all internal values of the matrix by taking the maximum value of
    # the top 'cell', the left 'cell', and the top-left 'cell' summed with the
    # computed identity for its corresponding loci
    for i in range(1, len(cluster_ref)+1):
        for j in range(1, len(cluster_alt)+1):
            
            if create_strings:
                comparison, identity = old_compare_loci(cluster_ref[i-1], cluster_alt[j-1], debug, verbose)
            else:
                comparison, identity, mismatch_zones = compare_loci(cluster_ref[i-1], cluster_alt[j-1], debug, verbose)
            
            if debug:
                print("comparison identity score = " + str(identity))
                print(f"top value = {dyn_prog_matrix[i-1][j]}; left value = {dyn_prog_matrix[i][j-1]}; diagonal value (identity) = {dyn_prog_matrix[i-1][j-1] + identity}")
            
            dyn_prog_matrix[i][j] = max(dyn_prog_matrix[i-1][j],
                                        dyn_prog_matrix[i][j-1],
                                        dyn_prog_matrix[i-1][j-1] + identity)
    if debug:        
        print("complete dynamic programmation matrix :")
        for line in dyn_prog_matrix:
            print(str(line))
        
    # retrieve best match alignment through backtracking
    
    i = len(cluster_ref)-1
    j = len(cluster_alt)-1
    results = []
    if debug:
        print(f"i = {i}, j = {j}")
    
    
    # while we did not yet reach the first column or line...
    while i>=0 and j>=0:
        if debug:
            print(f"i = {i}, j = {j}")
        
        if create_strings:
            comparison, identity = old_compare_loci(cluster_ref[i-1], cluster_alt[j-1], False, False)    
        else:
            comparison, identity, mismatch_zones = compare_loci(cluster_ref[i-1], cluster_alt[j-1], False, False)
        
        if debug:
            print(f"top value : {dyn_prog_matrix[i][j+1]}, \nleft value : {dyn_prog_matrix[i+1][j]}, \ndiagonal value : {dyn_prog_matrix[i][j] + identity}")
            print(f"max value : {max(dyn_prog_matrix[i][j+1], dyn_prog_matrix[i+1][j], dyn_prog_matrix[i][j] + identity)}")
        
        # if the maximum value is the diagonal value + computed identity, we
        # add both locus informations to the results dictionary and 'progress'
        # to the top-left cell
        if max(dyn_prog_matrix[i][j+1],
            dyn_prog_matrix[i+1][j],
            dyn_prog_matrix[i][j] + identity) == dyn_prog_matrix[i][j]+identity:
            if create_strings:
                comparison, identity = old_compare_loci(cluster_ref[i-1], cluster_alt[j-1], False, False)    
                if comparison == '_':
                    mismatch_zones = '_'
                else:
                    mismatch_zones = '?'            
            else:
                comparison, identity, mismatch_zones = compare_loci(cluster_ref[i-1], cluster_alt[j-1], False, False)
                
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
                            "reference mRNA number" : num_mRNAs_ref,
                            "alternative mRNA number" : num_mRNAs_alt})
            print(f"{results[-1]['cluster name']}\t\t{results[-1]['reference']}\t\t{results[-1]['alternative']}\t\t\t{results[-1]['mismatch/match']}\t\t\t\t{results[-1]['identity']}%")
            i -= 1
            j -= 1
        
        # if the maximum value is the top value, we add the reference locus 
        # informations to the results dictionary and 'progress' to the top cell
        elif max(dyn_prog_matrix[i][j+1],
            dyn_prog_matrix[i+1][j],
            dyn_prog_matrix[i][j] + identity) == dyn_prog_matrix[i][j+1]:
            if create_strings:
                comparison, identity = old_compare_loci(cluster_ref[i-1], cluster_alt[j], False, False)   
                mismatch_zones = '_'             
            else:
                comparison, identity, mismatch_zones = compare_loci(cluster_ref[i-1], cluster_alt[j], False, False)
                
            # we retrieve the number of mRNAs of the reference locus
            num_mRNAs_ref = len(cluster_ref[i-1].mRNAs.keys())
            num_mRNAs_alt = '_'
            
            results.append({"reference" : cluster_ref[i-1].name,
                            "reference start" : cluster_ref[i-1].start,
                            "reference end" : cluster_ref[i-1].end,
                            "alternative" : "_",
                            "alternative start" : "_",
                            "alternative end" : "_",
                            "mismatch/match" : comparison,
                            "identity" : identity,
                            "mismatch zones" : mismatch_zones,
                            "cluster name" : cluster_name,
                            "reference mRNA number" : num_mRNAs_ref,
                            "alternative mRNA number" : num_mRNAs_alt})
            print(f"{results[-1]['cluster name']}\t\t{results[-1]['reference']}\t\t{results[-1]['alternative']}\t\t\t{results[-1]['mismatch/match']}\t\t\t\t{results[-1]['identity']}%")
            i -= 1
           
        # if the maximum value is the left value, we add the alternative locus 
        # informations to the results dictionary and 'progress' to the left cell
        else:
            if create_strings:
                comparison, identity = old_compare_loci(cluster_ref[i], cluster_alt[j-1], False, False)     
                mismatch_zones = '_'           
            else:
                comparison, identity, mismatch_zones = compare_loci(cluster_ref[i], cluster_alt[j-1], False, False)
                
            # we retrieve the number of mRNAs of the alternative locus
            num_mRNAs_ref = '_'
            num_mRNAs_alt = len(cluster_alt[j-1].mRNAs.keys())
            
            results.append({"reference" : "_",
                            "reference start" : "_",
                            "reference end" : "_",
                            "alternative" : cluster_alt[j-1].name,
                            "alternative start" : cluster_alt[j-1].start,
                            "alternative end" : cluster_alt[j-1].end,
                            "mismatch/match" : comparison,
                            "identity" : identity,
                            "mismatch zones" : mismatch_zones,
                            "cluster name" : cluster_name,
                            "reference mRNA number" : num_mRNAs_ref,
                            "alternative mRNA number" : num_mRNAs_alt})
            print(f"{results[-1]['cluster name']}\t\t{results[-1]['reference']}\t\t{results[-1]['alternative']}\t\t\t{results[-1]['mismatch/match']}\t\t\t\t{results[-1]['identity']}%")
            j -= 1
            
    # if we reached the first line but did not reach the first column,
    # we add the rest of the alternative loci to the results
    while i==-1 and j!=-1:
    
        # we retrieve the number of mRNAs of the alternative locus
        num_mRNAs_ref = '_'
        num_mRNAs_alt = len(cluster_alt[j-1].mRNAs.keys())

        results.append({"reference" : "_",
                        "reference start" : "_",
                        "reference end" : "_",
                        "alternative" : cluster_alt[j-1].name,
                        "alternative start" : cluster_alt[j-1].start,
                        "alternative end" : cluster_alt[j-1].end,
                        "mismatch/match" : [],
                        "identity" : 0.0,
                        "mismatch zones" : "_",
                        "cluster name" : cluster_name,
                        "reference mRNA number" : num_mRNAs_ref,
                        "alternative mRNA number" : num_mRNAs_alt})
        print(f"{results[-1]['cluster name']}\t\t{results[-1]['reference']}\t\t{results[-1]['alternative']}\t\t\t{results[-1]['mismatch/match']}\t\t\t\t{results[-1]['identity']}%")
        j -= 1
        
    # if we reached the first column but did not reach the first line,
    # we add the rest of the reference loci to the results
    while i!=-1 and j==-1:
    
        # we retrieve the number of mRNAs of the reference locus
        num_mRNAs_ref = len(cluster_ref[i-1].mRNAs.keys())
        num_mRNAs_alt = '_'
        
        results.append({"reference" : cluster_ref[i-1].name,
                        "reference start" : cluster_ref[i-1].start,
                        "reference end" : cluster_ref[i-1].end,
                        "alternative" : "_",
                        "alternative start" : "_",
                        "alternative end" : "_",
                        "mismatch/match" : [],
                        "identity" : 0.0,
                        "mismatch zones" : "_",
                        "cluster name" : cluster_name,
                        "reference mRNA number" : num_mRNAs_ref,
                        "alternative mRNA number" : num_mRNAs_alt})
        print(f"{results[-1]['cluster name']}\t\t{results[-1]['reference']}\t\t{results[-1]['alternative']}\t\t\t{results[-1]['mismatch/match']}\t\t\t\t{results[-1]['identity']}%")
        i -= 1    
            
    return results


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
    
    results_file.write("Cluster name, Reference locus,Alternative locus,Comparison matches,Comparison mismatches,Identity score (%),Reference start,Reference end, Alternative start, Alternative end, non-correspondance zones, reference mRNA number, alternative mRNA number\n")
        
    # annotation origin of each locus in the results
    # (first value : both,  second value : reference,  third value : alternative)
    locus_initial_annot = [0,0,0]
    
    for cluster in results:
        for loc in cluster:
            # convert the mismatch zones so that commas don't modify 
            # the csv structure
            mismatch_zones = ""
            for zone in loc['mismatch zones']:
                zone += " // "
                mismatch_zones += zone
            mismatch_zones = mismatch_zones[:-4]
            
            # if no comparison was done for the loci, write '_' instead of 
            # the comparison values
            if loc['mismatch/match'] == []:
                results_file.write(f"{loc['cluster name']},{loc['reference']},{loc['alternative']},_,_,{loc['identity']},{loc['reference start']},{loc['reference end']},{loc['alternative start']},{loc['alternative end']}, _, {loc['reference mRNA number']}, {loc['alternative mRNA number']}\n")       
                if loc['reference'] == '_':
                    locus_initial_annot[2] += 1
                else:                       
                    locus_initial_annot[1] += 1
            else:
                results_file.write(f"{loc['cluster name']},{loc['reference']},{loc['alternative']},{loc['mismatch/match'][1]},{loc['mismatch/match'][0]},{loc['identity']},{loc['reference start']},{loc['reference end']},{loc['alternative start']},{loc['alternative end']}, {mismatch_zones}, {loc['reference mRNA number']}, {loc['alternative mRNA number']}\n")
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
# @return Returns a list of lists of dictionaries describing the 
# comparison of the structure identity between the loci of each annotation 
def annotation_comparison(ref_path, alt_path, debug=False, verbose=False, create_strings=False):

    # get all annotation files and generate the annotation data structure
    ref_annotations = get_gff_borders(ref_path, debug, verbose)
    alt_annotations = get_gff_borders(alt_path, debug, verbose)
    
    # get the order of the loci of both annotations
    locus_order = annotation_sort(ref_annotations, alt_annotations, debug, verbose)
    
    # construct clusters of overlapping loci
    clusters = construct_clusters(ref_annotations, alt_annotations, locus_order, debug, verbose)
    
    results = []
    print("\nCluster name\tReference_Locus\t\tAlternative_Locus\t\tComparison[mismatch, match]\t\tIdentity_Score\n")
    for cluster_id, cluster in clusters.clusters.items():
        results.append(annotation_match(cluster["ref"], cluster["alt"], cluster_id, create_strings, debug, verbose))
        
    write_results(results, debug, verbose)
    
    return results


def usage():
    
    # displayed when '-h' or '--help' is given, or when an invalid script
    # call happens
    print("Syntax : path/to/main.py [ -h/--help -d/--debug -v/--verbose -o/--old_version ] [ -r/--reference <reference_file_path> ] [ -a/--alternative <alternative_file_path> ] ")
    

def main():
    
    # we retrieve all script call options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hdvor:a:", ["help", "debug", "verbose", "old_version", "reference=", "alternative="])
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
        elif o in ("-r", "--reference"):
            ref_path = a
        elif o in ("-a", "--alternative"):
            alt_path = a
        else:
            assert False, "unhandled option"
            
    # call of the annotation_comparison function
    return annotation_comparison(ref_path, alt_path, debug, verbose, create_strings)
    
if __name__ == "__main__":
    main()
    

    
    