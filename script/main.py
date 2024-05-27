#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 11:57:29 2024

@author: vetea
"""

##@package main
# This script is used to compute the distance between two structural annotations of a same genome, one reference and one alternative annotation. 
# It expects as input the paths to the annotation files (in GFF format), 
# displays the computed distances between all annotation pairs, and returns a dictionary of lists of
# lists detailing the matchs/mismatchs between the two annotations' structure string


import getopt
import sys
import os
import pprint


#TODO documentation
class Clusters:
    def __init__(self):
        self.clusters = {}
        
    def clusters(self):
        return self.clusters
    
    def get_mRNAs(self, loc_id, ref):
        list_mRNAs = []
        for loc in self.clusters[loc_id][ref]:
            list_mRNAs.append(loc.mRNAs)
        return list_mRNAs
    
#TODO documentation
class Locus:
    def __init__(self, name="", mRNAs=None, start=-1, end=-1, direction=""):
        if not mRNAs:
            mRNAs = {}
        
        self.name = name
        self.mRNAs = mRNAs.copy()
        self.start = start
        self.end = end
        self.direction = direction
        
    def mRNAs(self):
        return self.mRNAs
    
    def contain_mrnas(self, **mrnas):
        for mrna_name, positions_list in mrnas.items():
            if mrna_name not in self.mRNAs:
                return False
            return self.mRNAs[mrna_name] == positions_list

    def show_init(self):
        return f"Locus(name='{self.name}', mRNAs={self.mRNAs}, start={self.start}, end={self.end}, direction='{self.direction}'"


## This function retrieves and returns the id of the structure described from a line read from a GFF file.
#
# @param line The line read from the file (string)
#
# @param debug If True, triggers display of many messages intended for debugging the program. Default is 'False'
#
# @param verbose If True, triggers display of more information messages. Default is 'False'
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
    
    # for each character from the first, we add it to the locus_id if it is not in a list of special characters. When the first special character is encountered, stop the loop and return the locus_id 
    while id_and_rest[i] not in [",", "?", ";", ":", "/", "!", "*", "$", "%", "+", "@", "#", "~", "&", "\n", "\t"] :
        
        structure_id += id_and_rest[i]
        i += 1
         
    if debug:
        print(f"Structure ID = {structure_id}")
        
    return structure_id


## This function retrieves and returns the parent id of the structure described from a line read from a GFF file.
#
# @param line The line read from the file (string)
#
# @param debug If True, triggers display of many messages intended for debugging the program. Default is 'False'
#
# @param verbose If True, triggers display of more information messages. Default is 'False'
#
# @remark This function expects the file to be in GFF format, and the parent id field to come after the id field
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
    
    # for each character from the first, we add it to the locus_id if it is not in a list of special characters. When the first special character is encountered, stop the loop and return the locus_id 
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
# @see Locus()
#
def get_gff_borders(path, debug=False, verbose=False):
    
    try:    
        file = open(path, "r") # the file to read
    except FileNotFoundError:
        print(f"\nget_gff_borders() function error : file '{path}' does not exist or cannot be read from\n")
        sys.exit(2)
    
    try:
        log = open("./results/log.txt", "w") # file in which to write structure errors
    except FileNotFoundError:
        os.mkdir("./results/") # create 'results' subdirectory
        log = open("./results/log.txt", "w") # file in which to write structure errors
        
    loci = {}
    locus_id = "" # the current gene being analyzed
    mRNA_id = "" # the current mRNA being analyzed
    line_index = 1 # number of the file line currently read
    locus = Locus()
    
    for line in file:
        
        # if we encounter a new gene, we get its ID and create a key in 'borders' with a basic list
        if str(line.split("\t")[2]) == "gene" and line.split("\t")[0] != "contig": 
            
            if locus_id != "" and locus.mRNAs[mRNA_id] == []: # if there was a previous gene, but its borders list is empty, return an error
                print(f"\nLine {line_index} = get_gff_borders() function error : no coding sequence (CDS) could be found for the previous mRNA '{mRNA_id}'\n")
                sys.exit(1)
            
            if locus.mRNAs != {}:
                loci[locus_id] = locus
            del locus
            locus = Locus()
            
            locus_id = get_structure_id(line, debug, verbose)
            locus.name = locus_id
            
            start = int(line.split("\t")[3])
            locus.start = start
            end = int(line.split("\t")[4])
            locus.end = end
            
            if locus.end < locus.start:
                locus.direction = "reverse"
            else:
                locus.direction = "direct"
            
            if verbose :
                print("\nReading the locus " + locus_id)

        # if we encounter a new gene, we get its ID and create a key in 'borders' with a basic list
        if str(line.split("\t")[2]) == "mRNA" and line.split("\t")[0] != "contig": 
            
            mRNA_id = get_structure_id(line, debug, verbose)
            locus.mRNAs[mRNA_id] = []
            if verbose :
                print("\nReading mRNA " + locus_id)

        # if we encounter a CDS, we add its start and end positions to the corresponding gene key in 'borders'
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
                
        line_index += 1
                
    if locus_id == "":
        print(f"\nget_gff_borders() function error : no locus could be found in file '{path}'\n")
        sys.exit(1)
        
    file.close()
    log.close()
    loci[locus_id] = locus
    
    # we return the entire dictionary with all loci
    return loci


## This function creates a list of all the locus coordinates for both annotations and sorts it in ascending order by their lower bound position. 
#
# @param dict_ref Dictionary containing all loci of the reference annotation, as returned by the 'get_gff_borders' function
#
# @param dict_alt Dictionary containing all loci of the alternative annotation, as returned by the 'get_gff_borders' function
#
# @param debug If True, triggers display of many messages intended for debugging the program. Default is 'False'
#
# @param verbose If True, triggers display of more information messages. Default is 'False'
#
# @see get_gff_borders()
#
# @return Returns a list of tuples containing the lower and upper bounds of each locus, its locus ID, and a boolean indicating if the locus was retrieved from the reference (True) or the alternative (False)
def annotation_sort(dict_ref, dict_alt, debug=False, verbose=False):
    
    if verbose:
        print("\nConstructing the locus order list of the two annotations")
 
    locus_order = []
    
    # get all reference loci bounds and locus ids as tuples in a list
    for locus_id, locus in dict_ref.items():
        locus_order.append((locus.start, locus.end, locus.name, True)) # 'True' indicates the tuple is from the reference
 		
    # get all alternative loci bounds and locus ids as tuples in a list
    for locus_id, locus in dict_alt.items():
        locus_order.append((locus.start, locus.end, locus.name, False)) # 'False' indicates the tuple is from the alternative
        
    if verbose:
        print("\nSorting the locus order list")
 	
    # sorts the list using the first value of each tuple (lower bound of the locus)
    locus_order.sort()
    
    if debug:
        print(f"\nlocus order list = {locus_order}")
     
    return locus_order


#TODO documentation
def construct_clusters(dict_ref, dict_alt, locus_order, debug=False, verbose=False):

    # initialisation of the dictionary keeping track of what loci have already been fused
    already_grouped = []
    
    clusters = Clusters()
    
    # for each locus in the loci list...
    for i in range(len(locus_order)-1):   
        locus_id = locus_order[i][2] # identifier of the current locus
        locus_borders = [locus_order[i][0], locus_order[i][1]] # lower and upper bounds of the current locus search
        locus_is_ref = locus_order[i][3]
        if debug:
            print(f"\nlocus borders = {locus_borders}")
            
        if locus_is_ref:
            if locus_id + "_ref" not in already_grouped:
                cluster_name = "cluster " + str(i) # name of the constructed cluster to add to the 'Clusters' class
                
                if cluster_name not in clusters.clusters:
                    clusters.clusters[cluster_name] = {"ref" : [],
                                                   "alt" : []} # initialize cluster in 'Clusters' class
                clusters.clusters[cluster_name]["ref"].append(dict_ref[locus_id])
        else:
            if locus_id + "_alt" not in already_grouped:
                cluster_name = "cluster " + str(i) # name of the constructed cluster to add to the 'Clusters' class
                
                if cluster_name not in clusters.clusters:
                    clusters.clusters[cluster_name] = {"ref" : [],
                                                   "alt" : []} # initialize cluster in 'Clusters' class
                clusters.clusters[cluster_name]["alt"].append(dict_alt[locus_id])
            
        # while we did not reach the end of the list and a 'stop signal' (j=-1) is not given, for each locus after the current one...
        j = 1
        while j != -1 and j <= len(locus_order)-i-1:
            next_locus_id = locus_order[i+j][2] # identifier of the current locus
            next_lower_border = locus_order[i+j][0] # lower bound of the locus pointed to by 'j'
            next_is_ref = locus_order[i+j][3] # boolean indicating if locus pointed by 'j' is from the reference
            if next_is_ref:
                is_ref_marker = "_ref"
            else:
                is_ref_marker = "_alt"
            
            if debug:
                print(f"\nj = {j}")
                print(f"next lower border = {next_lower_border}")
                print(f"already grouped : {already_grouped}")
            if next_locus_id + is_ref_marker not in already_grouped:
                if next_lower_border <= locus_borders[1]: # if the lower bound of the locus pointed by 'j' is inside the current locus' bounds and hasn't yet been included in a cluster...
                    if debug:
                        print("next locus found in current borders")
                    new_upper_border = locus_order[i+j][1]
                    locus_borders[1] = new_upper_border # upper bound of the current locus search is extended to upper bound of the 'j' locus               
                    if debug:
                        print(f"new locus borders = {locus_borders}")
                    
                    if next_is_ref:
                        clusters.clusters[cluster_name]["ref"].append(dict_ref[next_locus_id])
                        already_grouped.append(next_locus_id + is_ref_marker)
                    else:
                        clusters.clusters[cluster_name]["alt"].append(dict_alt[next_locus_id])
                        already_grouped.append(next_locus_id + is_ref_marker)
                    j += 1
                        
                else: # if no locus is found in the current locus search's bounds, a 'stop signal' is given to stop the 'j' loop
                    if debug:
                        print("next locus not found in current borders, continuing")
                    j = -1
            else:
                if debug:
                    print(f"{next_locus_id} already grouped, skipping")
                new_upper_border = locus_order[i+j][1]
                locus_borders[1] = new_upper_border # upper bound of the current locus search is extended to upper bound of the 'j' locus               
                if debug:
                    print(f"new locus borders = {locus_borders}")               
                j += 1
                         
        if debug:
            print(clusters.clusters[cluster_name]["ref"])
            print(clusters.clusters[cluster_name]["alt"])
    print("\n")
    
    return clusters


## This function retrieves all the CDS coordinates from the given lists of coordinates of the reference (@param ref) and of the alternative (@param alt) and includes them in a unique list of coordinates. The coordinates are sorted in ascending order.
#
# @param ref List of CDS coordinates of the reference annotation returned
#
# @param alt List of CDS coordinates of the alternative annotation returned
#
# @param debug If True, triggers display of many messages intended for debugging the program. Default is 'False'
#
# @param verbose If True, triggers display of more information messages. Default is 'False'
#
# @return Returns a list of coordinates compiling all coordinates from both initial lists in ascending order
def get_area_bounds(ref, alt, debug=False, verbose=False):
    
    if verbose:
        print("\nCreating comparison areas for the two annotations")
    
    bounds=[] # the return list
    i=0
    j=0    
    
    # while we did not reach the end of the coordinates lists...
    while i <= len(ref)-1 and j <= len(alt)-1:
        
        if debug:
            print(f"i = {i}, j = {j}")
            print(f"minimum of ref[i] ({ref[i]}) and alt[j] ({alt[j]}) = {min(ref[i], alt[j])}")
        
        # we add the next coordinate to the result list
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

## This function indicates for each couple of bounds in the given list of area bounds (@param area_bounds), if they delimit an area which includes a CDS from the given CDS coordinates list (@param cds_bounds). It is used during the comparison of areas in compare_loci() to know if the reference or alternative have a CDS in the area
#
# @see compare_loci()
#
# @param cds_bounds List of CDS coordinates for an annotation
#
# @param area_bounds List of area bounds
#
# @param debug If True, triggers display of many messages intended for debugging the program. Default is 'False'
#
# @param verbose If True, triggers display of more information messages. Default is 'False'
#
# @return Returns a list indicating for each couple of bounds if they include a CDS ('True') or not ('False')
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


#TODO change documentation
## This function compares two annotations' loci returned by the function get_gff_borders() and creates for each pair of reference-alternative mRNAs a comparison list detailing the identities and differences between the two annotations's codon position structure. It returns a tuple containing the comparison list giving the highest identity and the computed identity level.
#
# @see get_gff_borders()
#
# @param ref_locus The reference annotation's locus (Locus class instance)
#
# @param alt_locus The alternative annotation's locus (Locus class instance)
#
# @param debug If True, triggers display of many messages intended for debugging the program. Default is 'False'
#
# @param verbose If True, triggers display of more information messages. Default is 'False'
#
# @return Returns a tuple containign a list indicating the number of codon position mismatches (first position) and matches (second position) between the two border lists giving the maximum indentity between all mRNAs, and the identity level computed from this list
#
def compare_loci(ref_locus, alt_locus, debug=False, verbose=False):
    final_comparison = []
    final_identity = 0.0
    final_mismatch_zones = []
    
    for mRNA_ref_id, mRNA_ref in ref_locus.mRNAs.items():
        
        for mRNA_alt_id, mRNA_alt in alt_locus.mRNAs.items():
            if debug:
                print(f"mRNA_ref name = {mRNA_ref_id}")
                print(f"mRNA_alt_ name = {mRNA_alt_id}")
                print(f"mRNA_ref = {mRNA_ref}")
                print(f"mRNA_alt = {mRNA_alt}")
            
            # we retrieve the bounds of all areas delimited by all the CDS coordinates of both border lists
            area_bounds = get_area_bounds(mRNA_ref, mRNA_alt, debug, verbose)
            
            # we retrieve the list indicating the areas which include or not a CDS for both annotations
            
            if verbose:
                print("\nEvaluating presence of CDS in the comparison areas")
        
            ref_in_CDS = is_in_cds(mRNA_ref, area_bounds, debug, verbose)
            alt_in_CDS = is_in_cds(mRNA_alt, area_bounds, debug, verbose)
            if debug:
                print("Ref_in_CDS = " + str(ref_in_CDS))
                print("Alt_in_CDS = " + str(alt_in_CDS))
            
            codon_position_ref=0
            codon_position_alt=0
            comparison = [0,0] # return list. First value is mismatches, second value is identities
            mismatch_zones = []
            bound_id = 0
            prev_bounds = area_bounds[0]
            
            # for each comparison area delimited by the bounds in the list 'area_bounds'...
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
                    
                # if both annotations are in a CDS and have the same codon position, then all codon positions for the rest of the area will be identical, so we add the length of the area to the second value of the comparison list and update both codon positions
                if(ref_in_CDS[bound_id] and alt_in_CDS[bound_id] and codon_position_alt == codon_position_ref):
                    if debug:
                        print(f"Identical codon positions for the area, adding {bound-prev_bounds} to match values")
                    
                    comparison[1] += bound-prev_bounds
                    codon_position_alt = (codon_position_alt + (bound-prev_bounds))%3 
                    codon_position_ref = (codon_position_ref + (bound-prev_bounds))%3
                    
                # if both annotations are in a CDS but don't have the same codon position, then all codon positions for the rest of the area will be different, so we add the length of the are to the first value of the comparison list and update both codon positions
                elif(ref_in_CDS[bound_id] and alt_in_CDS[bound_id] and codon_position_alt != codon_position_ref):
                    if debug:
                        print(f"Different codon positions for the area, adding {bound-prev_bounds} to mismatch values")
                        
                    comparison[0] += bound-prev_bounds
                    codon_position_alt = (codon_position_alt + (bound-prev_bounds))%3 
                    codon_position_ref = (codon_position_ref + (bound-prev_bounds))%3 
                    mismatch_zones.append(f"{prev_bounds}-{bound}")
                
                # if only one annotation has a CDS in the comparison area, we add to the first value of the comparison list and update only one codon position
                
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
                        
                # the case of both annotations being outside of a CDS is not used in the computation of global loci identity, and is ignored (but mismatch zone borders are still collected)
                else:
                    if debug:
                        print(f"Reference and alternative are not in CDS for the area, ignoring {bound-prev_bounds} area")
                    
               
                prev_bounds = bound
                
                bound_id += 1
                
            # when the loci are on the reverse strand, matches and mismatches are counted in the negatives, so we convert them using the absolute value
            comparison[0] = abs(comparison[0])
            comparison[1] = abs(comparison[1])    
            
            if verbose:
                print(f"\nResult of the comparison of the locus : {comparison[1]} matches and {comparison[0]} mismatches\n")
                
            identity = comparison[1] / (comparison[0] + comparison[1])
        
            if identity > final_identity:
                final_comparison = comparison
                final_identity = identity
                final_mismatch_zones = mismatch_zones
    
    final_identity = round(final_identity * 100, 1)
    return (final_comparison, final_identity, final_mismatch_zones)


## This function expects a list of all CDS coordinates (start and end) of a locus. It returns a list indicating the start of the locus as first value and a string describing the codon position of each nucleotide (1,2,3, or 0 in the case of a non-CDS nucleotide) of the locus/gene as second value
#
# @param borders The list containing all start-end coordinates of the annotation's locus' CDS
#
# @see get_gff_borders()
#
# @return Returns a list describing the start of the locus and the annotation structure of the locus
def create_vectors(borders, debug=False, verbose=False):
    
    vector = [0, ""] # this variable takes in the strings of gene annotation structure for each gene
    
    if verbose:
        print("\nConverting the coordinates of the locus into a structure string")
    
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
        
        # we get the start (the last CDS value since the locus is reversed) of the locus
        vector[0] = borders[-1]
    
    # this variable takes the codon position of the next CDS nucleotide and loops
    # between the values 1, 2, and 3
    codon_pos = 1 
    
    # this variable indicates if we are in an exon/CDS or not.
    # it is incremented at each transition between annotations (each element in the 'borders' list)
    # to represent the exon-intron change along the gene
    in_exon = 1
    
    if borders[1] > borders[0]:
    
        # for each coordinate indexed for this locus in 'borders'
        for i in range(len(borders)-1):
        
            # if we are in a CDS, we append the numbers 1, 2, and 3 (with looping) 
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

        # for each coordinate indexed for this locus in 'borders' in reverse order
        for i in range(len(borders)-1, 0, -1):
        
            # if we are in a CDS, we append the numbers 1, 2, and 3 (with looping) 
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


## This function expects two structure strings corresponding to two annotations of the same genome, and compares each position of each string to return a list of mismatches (first value of the return list) and matches (secodn value of the return list)
#
# @param borders_vector_ref List of start position and vector of the locus of the first annotation
#
# @param borders_vector_alt List of start position and vector of the locus of the second annotation
#
# @see create_vectors()
#
# @return Returns a list describing the matchs/mismatchs for each position of the strings
#
# @remark This function doesn't expect any annotation to be a 'reference'
def old_compare_loci(ref_locus, alt_locus, debug=False, verbose=False): 
    final_comparison = []
    final_identity = 0.0
    
    for mRNA_ref_id, mRNA_ref in ref_locus.mRNAs.items():
        
        for mRNA_alt_id, mRNA_alt in alt_locus.mRNAs.items():
            vector_ref = create_vectors(mRNA_ref, debug, verbose)
            vector_alt = create_vectors(mRNA_alt, debug, verbose)
            comp_list = [0,0]
            if debug:
                print(f"vector_ref = {vector_ref}")
                print(f"vector_alt = {vector_alt}")
            
            # we get the minimum start positions, maximum end positions, and position difference of the loci
            minv = min(vector_ref[0], vector_alt[0]) # minimum start position
            diff = abs(vector_ref[0] - vector_alt[0]) # difference between the start positions
            if debug:
                print(f"minv = {minv}")
                print(f"diff = {diff}")
            
            # if the two loci don't start at the same position
            if vector_ref[0] != vector_alt[0]:
                if debug:
                    print("The two loci do not start at the same position")
                
                # if the locus 'a' starts before the start of locus 'b'
                if minv == vector_ref[0]:
                    if debug:
                        print("Locus 'a' start before locus 'b'")
                    
                    # for every comparison of the numbers at each position in the two strings, we increment by one the match or mismatch value of the return list. We account for the difference in start positions by adding the difference to the locus 'a' codon position retrieval
                    for i in range( min(len(vector_ref[1]), len(vector_alt[1])) + diff ):
                        
                        # try to know if the current position is in a CDS in the reference or alternative. If we are outside the string of an annotation the value 'False' is assigned
                        try:
                            ref_in_cds = vector_ref[1][i] in ("1", "2", "3") 
                        except IndexError:
                            ref_in_cds = False
                            
                        if i >= diff:
                            alt_in_cds = vector_alt[1][i-diff] in ("1", "2", "3")
                        else:
                            alt_in_cds = False
                        
                        if debug:
                            print("Range of loop = {min(len(vector_ref[1]), len(vector_alt[1])) + diff}")
                            print(f"i = {i}")
                            print(f"ref_in_cds = {ref_in_cds}")
                            print(f"alt_in_cds = {alt_in_cds}")
                        
                        # if we are outside the coordinates of locus 'b' and locus 'a' has a CDS in this position, we increment the mismatch count; else we do nothing
                        if i<diff and ref_in_cds:
                            comp_list[0] += 1
                        
                        # if we are outside the coordinates of locus 'a' and locus 'b' has a CDS in this position, we increment the mismatch count; else we do nothing
                        elif i >= len(vector_ref[1]) and alt_in_cds:
                            comp_list[0] += 1
                        
                        # if we are in both coordinates, we increment the match count if both have the same non-zero value at the current position, else we increment the mismatch value if they don't have both '0' as a value
                        else:
                            if ref_in_cds and vector_ref[1][i] == vector_alt[1][i-diff]:
                                comp_list[1] +=1
                            elif ref_in_cds or alt_in_cds:
                                comp_list[0] += 1
                
                # if the locus'a' starts after the start of locus 'b'
                elif minv == vector_alt[0]:
                    if debug:
                        print("Locus 'b' starts before locus 'a'")
                    
                    # for every comparison of the numbers at each position in the two strings, we increment by one the match or mismatch value of the return list. We account for the difference in start positions by adding the difference to the locus 'b' codon position retrieval
                    for i in range( min(len(vector_ref[1]), len(vector_alt[1])) + diff ):
                        
                        # try to know if the current position is in a CDS in the reference or alternative. If we are outside the string of an annotation the value 'False' is assigned
                        if i >= diff:
                            ref_in_cds = vector_ref[1][i-diff] in ("1", "2", "3")
                        else:
                            ref_in_cds = False
                            
                        try:
                            alt_in_cds = vector_alt[1][i] in ("1", "2", "3")
                        except IndexError:
                            alt_in_cds = False
                        
                        if debug:
                            print("Range of loop = {min(len(vector_ref[1]), len(vector_alt[1])) + diff")
                            print(f"i = {i}")
                            print(f"ref_in_cds = {ref_in_cds}")
                            print(f"alt_in_cds = {alt_in_cds}")
                            
                        # if we are outside the coordinates of locus 'a' and locus 'b' has a CDS in this position, we increment the mismatch count; else we do nothing
                        if i<diff and alt_in_cds:
                            comp_list[0] += 1
                            
                        # if we are outside the coordinates of locus 'b' and locus 'a' has a CDS in this position, we increment the mismatch count; else we do nothing
                        elif i >= len(vector_alt[1]) and ref_in_cds:
                            comp_list[0] += 1
                    
                        # if we are in both coordinates, we increment the match count if both have the same non-zero value at the current position, else we increment the mismatch value if they don't have both '0' as a value
                        else:
                            if ref_in_cds and vector_ref[1][i-diff] == vector_alt[1][i]:
                                comp_list[1] +=1
                            elif ref_in_cds or alt_in_cds:
                                comp_list[0] += 1
            
            # if the two loci start at the same position
            else:
                if debug:
                    print("The two loci start at the same position")
            
                # for every comparison of the numbers at each position in the two strings, we increment by the match value of the return list of the two string values at this position are equal, or the mismatch value if they are not or if only one annotation has a CDS at this position
                for i in range(len(vector_ref[1])):
                    
                    if i < len(vector_ref[1]):
                        ref_in_cds = vector_ref[1][i-diff] in ("1", "2", "3")
                    else:
                        ref_in_cds = False
                        
                    if i < len(vector_alt[1]):
                        alt_in_cds = vector_alt[1][i-diff] in ("1", "2", "3")
                    else:
                        alt_in_cds = False
                
                    if debug:
                        print(f"\nRange of loop = {len(vector_ref[1])}")
                        print(f"i = {i}")
                        print(f"ref_in_cds = {ref_in_cds}")
                        print(f"alt_in_cds = {alt_in_cds}")
                    
                    if ref_in_cds and alt_in_cds and vector_ref[1][i] == vector_alt[1][i]:
                        comp_list[1] += 1
                    elif ref_in_cds or alt_in_cds:
                        comp_list[0] += 1
        
            identity = comp_list[1] / (comp_list[0] + comp_list[1])
        
            if identity > final_identity:
                final_comparison = comp_list
                final_identity = identity
    
    final_identity = round(final_identity * 100, 1)
        
    if verbose:
        print(f"\nResult of the comparison of the locus : {comp_list[1]} matches and {comp_list[0]} mismatches" )
    
    return (final_comparison, final_identity)
    

#TODO documentation
#TODO add parameter to select old or new version of compare_loci()
def annotation_match(cluster_ref, cluster_alt, create_strings, debug=False, verbose=False):
    
    if verbose:
        print("\nmatching annotations loci with each other")
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
            
    # compute all internal values of the matrix
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
    
    print("\nReference_Locus\t\tAlternative_Locus\t\tComparison[mismatch, match]\t\tIdentity_Score\n")
    
    while i>=0 and j>=0:

        if create_strings:
            comparison, identity = old_compare_loci(cluster_ref[i-1], cluster_alt[j-1], False, False)    
        else:
            comparison, identity, mismatch_zones = compare_loci(cluster_ref[i-1], cluster_alt[j-1], False, False)
        
        if debug:
            print(f"top value : {dyn_prog_matrix[i][j+1]}, \nleft value : {dyn_prog_matrix[i+1][j]}, \ndiagonal value : {dyn_prog_matrix[i][j] + identity}")
            print(f"max value : {max(dyn_prog_matrix[i][j+1], dyn_prog_matrix[i+1][j], dyn_prog_matrix[i][j] + identity)}")
        
        if max(dyn_prog_matrix[i][j+1],
            dyn_prog_matrix[i+1][j],
            dyn_prog_matrix[i][j] + identity) == dyn_prog_matrix[i][j]+identity:
            if create_strings:
                comparison, identity = old_compare_loci(cluster_ref[i-1], cluster_alt[j-1], False, False)                
            else:
                comparison, identity, mismatch_zones = compare_loci(cluster_ref[i-1], cluster_alt[j-1], False, False)
            results.append({"reference" : cluster_ref[i-1].name,
                            "reference start" : cluster_ref[i-1].start,
                            "reference end" : cluster_ref[i-1].end,
                            "alternative" : cluster_alt[j-1].name,
                            "alternative start" : cluster_alt[j-1].start,
                            "alternative end" : cluster_alt[j-1].end,
                            "mismatch/match" : comparison,
                            "identity" : identity,
                            "mismatch zones" : mismatch_zones})
            print(f"{results[-1]['reference']}\t\t{results[-1]['alternative']}\t\t\t{results[-1]['mismatch/match']}\t\t\t\t{results[-1]['identity']}%")
            i -= 1
            j -= 1
        
        elif max(dyn_prog_matrix[i][j+1],
            dyn_prog_matrix[i+1][j],
            dyn_prog_matrix[i][j] + identity) == dyn_prog_matrix[i][j+1]:
            if create_strings:
                comparison, identity = old_compare_loci(cluster_ref[i-1], cluster_alt[j], False, False)                
            else:
                comparison, identity, mismatch_zones = compare_loci(cluster_ref[i-1], cluster_alt[j], False, False)
            results.append({"reference" : cluster_ref[i-1].name,
                            "reference start" : cluster_ref[i-1].start,
                            "reference end" : cluster_ref[i-1].end,
                            "alternative" : "_",
                            "alternative start" : "_",
                            "alternative end" : "_",
                            "mismatch/match" : comparison,
                            "identity" : identity,
                            "mismatch zones" : mismatch_zones})
            print(f"{results[-1]['reference']}\t\t{results[-1]['alternative']}\t\t\t{results[-1]['mismatch/match']}\t\t\t\t{results[-1]['identity']}%")
            i -= 1
            
        else:
            if create_strings:
                comparison, identity = old_compare_loci(cluster_ref[i], cluster_alt[j-1], False, False)                
            else:
                comparison, identity, mismatch_zones = compare_loci(cluster_ref[i], cluster_alt[j-1], False, False)
            results.append({"reference" : "_",
                            "reference start" : "_",
                            "reference end" : "_",
                            "alternative" : cluster_alt[j-1].name,
                            "alternative start" : cluster_alt[j-1].start,
                            "alternative end" : cluster_alt[j-1].end,
                            "mismatch/match" : comparison,
                            "identity" : identity,
                            "mismatch zones" : mismatch_zones})
            print(f"{results[-1]['reference']}\t\t{results[-1]['alternative']}\t\t\t{results[-1]['mismatch/match']}\t\t\t\t{results[-1]['identity']}%")
            j -= 1
            
        while i==-1 and j!=-1:
            if create_strings:
                comparison, identity = old_compare_loci(cluster_ref[i], cluster_alt[j-1], False, False)                
            else:
                comparison, identity, mismatch_zones = compare_loci(cluster_ref[i], cluster_alt[j-1], False, False)
            results.append({"reference" : "_",
                            "reference start" : "_",
                            "reference end" : "_",
                            "alternative" : cluster_alt[j-1].name,
                            "alternative start" : cluster_alt[j-1].start,
                            "alternative end" : cluster_alt[j-1].end,
                            "mismatch/match" : comparison,
                            "identity" : identity,
                            "mismatch zones" : mismatch_zones})
            print(f"{results[-1]['reference']}\t\t{results[-1]['alternative']}\t\t\t{results[-1]['mismatch/match']}\t\t\t\t{results[-1]['identity']}%")
            j -= 1
            
        while i!=-1 and j==-1:
            if create_strings:
                comparison, identity = old_compare_loci(cluster_ref[i-1], cluster_alt[j], False, False)                
            else:
                comparison, identity, mismatch_zones = compare_loci(cluster_ref[i-1], cluster_alt[j], False, False)
            results.append({"reference" : cluster_ref[i-1].name,
                            "reference start" : cluster_ref[i-1].start,
                            "reference end" : cluster_ref[i-1].end,
                            "alternative" : "_",
                            "alternative start" : "_",
                            "alternative end" : "_",
                            "mismatch/match" : comparison,
                            "identity" : identity,
                            "mismatch zones" : mismatch_zones})
            print(f"{results[-1]['reference']}\t\t{results[-1]['alternative']}\t\t\t{results[-1]['mismatch/match']}\t\t\t\t{results[-1]['identity']}%")
            i -= 1    
            
    return results


## This function writes the results of the annotation comparison retrieved from the identities dictionary to a new 'results.csv' file
#
# @param results A list of list of dictionaries containing results of the annotation comparison, as returned by annotation_match
#
# @param debug If True, triggers display of many messages intended for debugging the program. Default is 'False'
#
# @param verbose If True, triggers display of more information messages. Default is 'False'
#
# @see annotation_match()
#
def write_results(results, debug=False, verbose=False):
    
    try:
        results_file = open("./results/results.csv", "w") # file in which to write results
    except FileNotFoundError:
        os.mkdir("./results/") # create 'results' subdirectory
        results_file = open("./results/results.csv", "w") # file in which to write results
    
    results_file.write("Reference locus,Alternative locus,Comparison matches,Comparison mismatches,Identity score (%),Reference start,Reference end, Alternative start, Alternative end, non-correspondance zones\n")
        
    for cluster in results:
        for loc in cluster:
            # converting the mismatch zones so that commas don't modifiy 
            # the csv structure
            mismatch_zones = ""
            for zone in loc['mismatch zones']:
                zone += " // "
                mismatch_zones += zone
            mismatch_zones = mismatch_zones[:-4]
            
            
            if loc['mismatch/match'] == []:
                results_file.write(f"{loc['reference']},{loc['alternative']},_,_,{loc['identity']},{loc['reference start']},{loc['reference end']},{loc['alternative start']},{loc['alternative end']}, _\n")                                
            else:
                results_file.write(f"{loc['reference']},{loc['alternative']},{loc['mismatch/match'][1]},{loc['mismatch/match'][0]},{loc['identity']},{loc['reference start']},{loc['reference end']},{loc['alternative start']},{loc['alternative end']}, {mismatch_zones}\n")
        
    results_file.close()
    
    
## Main function of this program. Given a reference and alternative path, gets the corresponding GFF files and compares the two annotations to return their structure's identity level
#
# @param ref_path Path of the GFF file describing the reference annotation
#
# @param alt_path Path of the GFF file describing the aternative annotation
#
# @param debug If True, triggers display of many messages intended for debugging the program. Default is 'False'
#
# @param verbose If True, triggers display of more information messages. Default is 'False'
#
# @return Returns a dictionary of dictionaries of floats corresponding to the structure
# string identity between each locus of each annotation compared to those of the reference
#
# @remark Loci found in one annotation but not the other are ignored
def annotation_comparison(ref_path, alt_path, debug=False, verbose=False, create_strings=False):

    # get all annotation files and generate the annotation data structure
    ref_annotations = get_gff_borders(ref_path, debug, verbose)
    alt_annotations = get_gff_borders(alt_path, debug, verbose)
    
    # get the order of the loci of both annotations
    locus_order = annotation_sort(ref_annotations, alt_annotations, debug, verbose)
    
    # construct clusters of overlapping loci
    clusters = construct_clusters(ref_annotations, alt_annotations, locus_order, debug, verbose)
    
    results = []
    for cluster_id, cluster in clusters.clusters.items():
        results.append(annotation_match(cluster["ref"], cluster["alt"], create_strings, debug, verbose))
        
    write_results(results, debug, verbose)
    
    # identities = {}
    
    # # for each locus of the reference annotation...
    # for locus_id, locus in ref_annotations.items():
    
    #     # if the locus exists in the alternative annotation and if the two loci are on the same strand (same strand number)...
    #     if locus_id in alt_annotations and ref_annotations[locus_id].direction == alt_annotations[locus_id].direction:
        
    #         # construct the comparison matrix for the reference and alternative locus annotation...
              
    #         # using the new version
    #         if create_strings == False:
    #             if verbose:
    #                 print("\nStarting the comparison of the locus " + locus + " using the new program version\n")
    #             comp_res = compare_loci(ref_annotations[locus][1], alt_annotations[locus][1], debug, verbose)
                
    #         # using the old version
    #         else:
    #             if verbose:
    #                 print("\nStarting the comparison of the locus " + locus + " using the old program version (structure strings creation)\n")
    #             comp_res = old_compare_loci(ref_annotations[locus][1], alt_annotations[locus][1], debug, verbose)
            
    #         # compute identity from the matrix and index it in the identities dictionary
    #         identities[locus] = round( float(comp_res[1]) / (float(comp_res[1]) + float(comp_res[0]) ) * 100 , 1 )
            
    #         if verbose:
    #             print("Finished computing the identity of the annotations of locus " + locus + "\n")
             
    #     # if the two loci are on different strands or if the locus doesn't exist in the alternative annotation, we assign 0% identity to the locus
    #     else:
        
    #         identities[locus] = 0.0
            
    #         if verbose:
                
    #             if locus in alt_annotations:
    #                 print("\nAlternative locus is predicted on a different strand from the reference. Identity : 0%\n")
    #             else:
    #                 print("\nLocus is not predicted in alternative annotation. Identity: 0%\n")
                    
    # # for each locus of the alternative annotation...
    # for locus in alt_annotations:
        
    #     # if the locus doesn't exist in the reference, assign 0% identity
    #     if locus not in ref_annotations:
            
    #         if verbose:
    #             print(f"\nLocus {locus} found in alternative annotation, but not in reference. Identity: 0%\n")
    #         identities[locus] = 0.0
                
    # write_results(identities, debug, verbose)        
    
    # return identities


def usage():
    
    # displayed when '-h' or '--help' is given, or when an invalid script call happens
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
    # debug: display lots of informations on internal function variables and mecanisms, not intended to be used by the end user
    # verbose: display messages indicating which step the program is currently on, intended to be used when the program is called directly (not integrated in a pipeline or workflow)
    debug = False
    verbose = False
    
    # boolean indicating which version of the program to use: False = new version without any structure string creation, True = 'old' version with creation of structure strings to compare the loci of the annotations
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
    annotation_comparison(ref_path, alt_path, debug, verbose, create_strings)
    
    # print("\nResult of the comparison of all loci of both annotations :\n")
    # pprint.pprint(comparison)
    # print("\n\n")
    
if __name__ == "__main__":
    main()
    

    
    