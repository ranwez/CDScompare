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
import pprint


## This function expects a string corresponding to the file path of the GFF file to read, and returns
# a dictionary of lists with the keys corresponding to a locus identifier, and the values 
# corresponding to a list of a number indicating the strand supporting the locus (1 for direct 
# strand, 0 for reverse strand) and a list of the start and end position of each coding 
# sequence ('CDS') of the locus 
#
# @param path Path of the file to read
#
# @param debug If True, triggers display of many messages intended for debugging the program
#
# @param verbose If True, triggers display of more information messages. Default is 'False'
#
# @return Returns a dictionary of lists containing the number associated with the strand of the locus (1 : direct, 0: reverse) and a list containing the start and end coordinates of all CDS for each locus (identified with the locus ID)
def get_gff_borders(path, debug=False, verbose=False):
    
    try:    
        file = open(path, "r") # the file to read
    except FileNotFoundError:
        print(f"\nget_gff_borders() function error : file '{path}' does not exist or cannot be read from\n")
        sys.exit(2)
    
    borders = {} # this variable takes in the borders of each CDS of each gene
    locus_id = "" # the current gene being analyzed
    
    for l in file:

        # if we encounter a new gene, we get its ID and create a key in 'borders' with a basic list
        if str(l.split("\t")[2]) == "gene": 
            
            if locus_id != "":
                if borders[locus_id] == [1, []]:
                    print(f"\nget_gff_borders() function error : no coding sequence (CDS) could be found for the locus '{locus_id}'\n")
                    sys.exit(1)
            
            locus_id = l.split("\t")[8].split("\n")[0][3:]
            
            # the locus list is intialised with a strand number of '1', which is the corrected if necessary
            borders[locus_id] = [1, []]
            
            if verbose :
                print("\nReading the locus " + locus_id)

        # if we encounter a CDS, we add its start and end positions to corresponding gene key in 'borders'
        if str(l.split("\t")[2]) == "CDS":
            
            borders[locus_id][1].append(int(l.split("\t")[3]))
            borders[locus_id][1].append(int(l.split("\t")[4]))
            
            if verbose:
                print("Adding borders to " + locus_id + " : " + l.split("\t")[3] + ", " + l.split("\t")[4])
                
            # if the retrieved borders indicate the locus is on the reverse strand, we change the strand number to '0'
            if borders[locus_id][0] == 1 and borders[locus_id][1][1] < borders[locus_id][1][0]:
                
                borders[locus_id][0] = 0
                
    if borders == {}:
        print(f"\nget_gff_borders() function error : no locus could be found in file '{path}'\n")
        sys.exit(1)
        
    file.close()
    
    # we return the entire dictionary with all borders
    return borders


## This function creates a list of all the locus coordinates for both annotations and sorts it in ascending order by their lower bound position. 
#
# @param dict_ref Dictionary containing all loci of the reference annotation, as returned by the 'get_gff_borders' function
#
# @param dict_alt Dictionary containing all loci of the alternative annotation, as returned by the 'get_gff_borders' function
#
# @param debug If True, triggers display of many messages intended for debugging the program
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
    
    # get all reference loci bounds and locus ids
    for locus_id, borders in dict_ref.items():
        locus_order.append((borders[1][0], borders[1][-1], locus_id, True)) # 'True' indicates the tuple is from the reference
 		
    # get all alternative loci bounds and locus ids
    for locus_id, borders in dict_alt.items():
        locus_order.append((borders[1][0], borders[1][-1], locus_id, False)) # 'False' indicates the tuple is from the alternative
        
    if verbose:
        print("\nSorting the locus order list")
 	
    # sorts the list using the first value of each tuple (lower bound of the locus)
    locus_order.sort()
    
    if debug:
        print(f"\nlocus order list = {locus_order}")
     
    return locus_order


## This function merges two loci from an annotation dictionary by appending each element (coordinate) of the second locus to the list of coordinates of the first locus and deleting the locus fro mthe annotation
#
# @param annotation_dict The dictionary of the annotation to modifiy, as returned by the function 'get_gff_borders'
#
# @param loc_a_id Locus id of the locus into which the other is merged
#
# @param loc_b_id Locus id of the locus to merge into the other and delete
#
# @param debug If True, triggers display of many messages intended for debugging the program
#
# @param verbose If True, triggers display of more information messages. Default is 'False'
#
# @see get_gff_borders()
#
# @remark This function returns nothing and is expected to be used without affectation. This function modifies the original dictionary passed and deletes information from it.
#
# @return Returns nothing
def locus_append_delete(annotation_dict, loc_a_id, loc_b_id, debug=False, verbose=False):
    
    if verbose:
        print(f"\nFusing the locus {loc_b_id} into the locus {loc_a_id}")
    
    # append all coordinates of 'locus_b' to the list of 'locus a'
    for i in range(len(annotation_dict[loc_b_id][1])):
        
        annotation_dict[loc_a_id][1].append(annotation_dict[loc_b_id][1][i])
        
    # delete 'locus b' from the dictionary
    del annotation_dict[loc_b_id]


## This function fuses all overlapping loci of both annotation dictionaries by appending their coordinates.
#
# @param dict_ref Dictionary of the reference annotation, as returned by the function 'get_gff_borders'
#
# @param dict_alt Dictionary of the alternative annotation, as returned by the function 'get_gff_borders'
#
# @param locus_order List of tuples of lower and upper bounds, locus id, and boolean of reference origin for each locus of both annotations, as returned by the function 'annotation_sort'
#
# @see get_gff_borders()
# @see annotation_sort()
#
# @remark This function returns nothing and is expected to be used without affectation. This function modifies the original dictionaries passed and deletes information from them.
def fuse_superloci(dict_ref, dict_alt, locus_order, debug=False, verbose=False):

    # initialisation of the dictionary keeping track of what loci have already been fused
    already_fused = {'ref': {},
                     'alt': {}}

    # for each locus in the loci list...
    for i in range(len(locus_order)-1):   
        locus_id = locus_order[i][2] # identifier of the current locus
        is_ref = locus_order[i][3] # boolean indicating if current locus is from the reference
        locus_borders = [locus_order[i][0], locus_order[i][1]] # lower and upper bounds of the current locus
        if debug:
            print(f"\nlocus borders = {locus_borders}")
	 
        # while we did not reach the end of the list and a 'stop signal' (j=-1) is not given, for each locus after the current one...
        j = 1
        while j != -1 and j <= len(locus_order)-i-1:
            next_lower_border = locus_order[i+j][0] # lower bound of the locus pointed to by 'j'
            next_is_ref = locus_order[i+j][3] # boolean indicating if locus pointed by 'j' is from the reference
            if debug:
                print(f"\nj = {j}")
                print(f"next lower border = {next_lower_border}")
            
            if is_ref:
                if debug:
                    print("current is in ref")
                
                if next_lower_border <= locus_borders[1]: # if the lower bound of the locus pointed by 'j' is inside the current locus' bounds...
                    if debug:
                        print("next locus found in current borders")
                    new_upper_border = locus_order[i+j][1]
                    locus_borders[1] = new_upper_border # upper bound of the current locus is extended to upper bound of the 'j' locus               
                    
                    # if the locus pinted by 'j' is in the reference, it is merged with the current one, else we do nothing
                    if next_is_ref:
                        next_locus_id = locus_order[i+j][2]
                        
                        if next_locus_id not in already_fused['ref']:
                            if debug:
                                print("next is in ref, fusing with current locus")
                            locus_append_delete(dict_ref, locus_id, next_locus_id, debug, verbose)
                            already_fused['ref'][next_locus_id] = True
                        else:
                            if debug:
                                print("next was already fused, continuing...")
                    else:
                        if debug:
                            print("next is not in ref")
                    if debug:
                        print(f"new locus borders = {locus_borders}")
                    j += 1
                        
                else: # if no locus is found in the current locus bounds, a 'stop signal' is given to stop the 'j' loop
                    if debug:
                        print("next locus not found in current borders, continuing")
                    j = -1
                
            else:
                if debug:
                    print("current is in alt")
                
                if next_lower_border <= locus_borders[1]: # if the lower bound of the locus pointed by 'j' is inside the current locus' bounds...
                    if debug:
                        print("next locus found in current borders")
                    new_upper_border = locus_order[i+j][1]
                    locus_borders[1] = new_upper_border # upper bound of the current locus is extended to upper bound of the 'j' locus
                    
                    # if the locus pinted by 'j' is in the alternative, it is merged with the current one, else we do nothing
                    if not next_is_ref:
                        next_locus_id = locus_order[i+j][2]
                        
                        if next_locus_id not in already_fused['alt']:
                            if debug:
                                print("next is in alt, fusing with current locus")
                            locus_append_delete(dict_alt, locus_id, next_locus_id, debug, verbose)
                            already_fused['alt'][next_locus_id] = True
                        else:
                            if debug:
                                print("next was already fused, continuing...")
                    else:
                        if debug:
                            print("next is not in alt")
                    if debug:
                        print(f"new locus borders = {locus_borders}")
                    j += 1
                else: # if no locus is found in the current locus bounds, a 'stop signal' is given to stop the 'j' loop
                    if debug:
                        print("next locus not found in current borders, continuing")
                    j = -1
                     
        if debug:
            print(dict_ref)
            print(dict_alt)
    print("\n")


## This function retrieves all the CDS coordinates from the given lists of coordinates of the reference (@param ref) and of the alternative (@param alt) and includes them in a unique list of coordinates. The coordinates are sorted in ascending order.
#
# @param ref List of CDS coordinates of the reference annotation returned
#
# @param alt List of CDS coordinates of the alternative annotation returned
#
# @param debug If True, triggers display of many messages intended for debugging the program
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
            
            bounds.append(alt[j])
        
    if j > len(alt)-1:
        
        for k in range(i, len(ref)):
            
            bounds.append(ref[i])
            
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
# @param debug If True, triggers display of many messages intended for debugging the program
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
        
        current_cds=[cds_bounds[cds_id],cds_bounds[cds_id+1]]
        lb= area_bounds[i]
        ub= area_bounds[i+1]
        
        if debug:
            print(f"\ni = {i}")
            print(f"cds_id = {cds_id}")
            print(f"current_cds = {current_cds}")
            print(f"Lower bound (lb) = {lb}")
            print(f"Upper bound (ub) = {ub}")

        # if the current CDS is included in the area bounds, 'True' is used to index it in the return list, else 'False'
        if(current_cds[0]<=lb and ub<=current_cds[1]):
            in_cds.append(True)
        else:
            in_cds.append(False)
            
        if debug:
            print(f"in_cds = {in_cds}")
            
        # if the current area reaches the end of the current CDS, we skip to the next CDS
        if(ub==current_cds[1]):
            cds_id+=2
            
        i+=1
        
    # if we reached the end of the CDS coordinates list, we continue until all area bounds are exhausted by appending 'False'
    while(i < len(area_bounds)-1):
        
        if debug:
            print(f"i (after end of CDS list) = {i}")
        
        in_cds.append(False)
        
        i += 1
        
    if debug:
        print(f"\nFinal in_cds = {in_cds}\n")
        
    return in_cds

## This function compares two annotations' loci from their CDS border list returned by the function get_gff_borders() and creates a comparison list detailing the identities and differences between the two annotations's codon position structure
#
# @see get_gff_borders()
#
# @param ref The reference annotation's border list
#
# @param alt The alternative annotation's border list
#
# @param debug If True, triggers display of many messages intended for debugging the program
#
# @param verbose If True, triggers display of more information messages. Default is 'False'
#
# @return Returns a list indicating the number of codon position mismatches (first position) and matches (second position) between the two border lists
#
def compare_loci(ref, alt, debug=False, verbose=False):
    
    # we retrieve the bounds of all areas delimited by all the CDS coordinates of both border lists
    area_bounds = get_area_bounds(ref, alt, debug, verbose)
    
    # we retrieve the list indicating the areas which include or not a CDS for both annotations
    
    if verbose:
        print("\nEvaluating presence of CDS in the comparison areas")

    ref_in_CDS = is_in_cds(ref, area_bounds, debug, verbose)
    alt_in_CDS = is_in_cds(alt, area_bounds, debug, verbose)
    
    if debug:
        print("Ref_in_CDS = " + str(ref_in_CDS))
        print("Alt_in_CDS = " + str(alt_in_CDS))
    
    codon_position_ref=0
    codon_position_alt=0
    
    comparison = [0,0] # return list. First value is mismatches, second value is identities
    bound_id=0
    prev_bounds=area_bounds[0]
    
    # for each comparison area delimited by the bounds in the list 'area_bounds'...
    for bound in area_bounds[1:] :
        
        if(prev_bounds>0):
                
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
            if(ref_in_CDS[bound_id] and alt_in_CDS[bound_id] and codon_position_alt==codon_position_ref):
                
                if debug:
                    print(f"Identical codon positions for the area, adding {bound-prev_bounds} to match values")
                
                comparison[1]+=bound-prev_bounds
                codon_position_alt = (codon_position_alt + (bound-prev_bounds))%3 
                codon_position_ref = (codon_position_ref + (bound-prev_bounds))%3
                
            # if both annotations are in a CDS but don't have the same codon position, then all codon positions for the rest of the area will be different, so we add the length of the are to the first value of the comparison list and update both codon positions
            elif(ref_in_CDS[bound_id] and alt_in_CDS[bound_id] and codon_position_alt!=codon_position_ref):
                
                if debug:
                    print(f"Different codon positions for the area, adding {bound-prev_bounds} to mismatch values")
                    
                comparison[0]+=bound-prev_bounds
                codon_position_alt = (codon_position_alt + (bound-prev_bounds))%3 
                codon_position_ref = (codon_position_ref + (bound-prev_bounds))%3 
            
            # if only one annotation has a CDS in the comparison area, we add to the first value of the comparison list and update only one codon position
            
            elif(ref_in_CDS[bound_id]):
                
                if debug:
                    print(f"Alternative is not in CDS for the area, adding {bound-prev_bounds} to mismatch values")
                    
                comparison[0]+=bound-prev_bounds
                codon_position_ref = (codon_position_ref + (bound-prev_bounds))%3 
                
            
            elif(alt_in_CDS[bound_id]):
                
                if debug:
                    print(f"Reference is not in CDS for the area, adding {bound-prev_bounds} to mismatch values")
                    
                comparison[0]+=bound-prev_bounds
                codon_position_alt = (codon_position_alt + (bound-prev_bounds))%3 
                
            # the case of both annotations being outside of a CDS is not used in the computation of global loci identity, and is ignored
       
        prev_bounds=bound
        
        bound_id += 1
        
    # when the loci are on the reverse strand, matches and mismatches are counted in the negatives, so we convert them using the absolute value
    comparison[0] = abs(comparison[0])
    comparison[1] = abs(comparison[1])    
    
    if verbose:
        print(f"\nResult of the comparison of the locus : {comparison[1]} matches and {comparison[0]} mismatches\n")
        
    return comparison


## This function expects a list of all CDS coordinates (start and end) of a locus. It returns a list indicating the start of the locus as first value and a string describing the codon position of each nucleotide (1,2,3, or 0 in the case of a non-CDS nucleotide) of the locus/gene as second value
#
# @param borders The list containing all start-end coordinates of the annotation's locus' CDS
#
# @see get_gff_borders()
#
# @return Returns a list describing the start of the locus and the annotation structure of the locus
def create_vectors(borders, debug=False, verbose=False):
    
    vectors = [0, ""] # this variable takes in the strings of gene annotation structure for each gene
    
    if verbose:
        print("\nConverting the coordinates of the locus into a structure string")
    
    # if the locus is on the direct strand
    if borders[1] > borders[0]:
        
        if debug:
            print(f"\nLocus is on direct strand, retrieving start of locus : {borders[0]}")
        
        # we get the start of the locus
        vectors[0] = borders[0]
        
    # if the locus is on the reverse strand
    else:
        
        if debug:
            print(f"\nLocus is on reverse strand, retrieving start of locus from end of coordinates list : {borders[-1]}")
        
        # we get the start (the last CDS value since the locus is reversed) of the locus
        vectors[0] = borders[-1]
    
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
                    vectors[1] += str(codon_pos)
                    
                    # we increment the codon position with looping
                    if codon_pos == 3:
                        codon_pos = 1
                    else:
                        codon_pos += 1
                        
                if debug:
                    print(f"New codon position = {codon_pos}")
                    print(f"New structure string : {vectors[1]}")
                
                
            # if we are not in a CDS, we append 0
            if in_exon % 2 == 0:
                
                if debug:
                    print(f"i = {i}")
                    print("in_exon = False (adding 0 to structure string)")
                
                # for each nucleotide between this coordinate and the next...
                for j in range( borders[i+1] - borders[i] ):                
                
                    vectors[1] += "0"
                    
                if debug:
                    print(f"New structure string : {vectors[1]}")
            
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
                    vectors[1] += str(codon_pos)
                    
                    # we increment the codon position with looping
                    if codon_pos == 3:
                        codon_pos = 1
                    else:
                        codon_pos += 1
                        
                if debug:
                    print(f"New codon position = {codon_pos}")
                    print(f"New structure string : {vectors[1]}")
                
            # if we are not in a CDS, we append 0
            if in_exon % 2 == 0:

                if debug:
                    print(f"i = {i}")
                    print("in_exon = False (adding 0 to structure string)")
                
                # for each nucleotide between this coordinate and the next...
                for j in range( borders[i-1] - borders[i] ):
                
                    vectors[1] += "0"
                    
                if debug:
                    print(f"New structure string : {vectors[1]}")
            
            in_exon += 1   
    
    if verbose:
        print("\nStructure string of the locus :\n" + vectors[1] + "\n")
        
    return vectors


## This function expects two structure strings corresponding to two annotations of the same genome, and compares each position of each string to return a list of mismatches (first value of the return list) and matches (secodn value of the return list)
#
# @param borders_loc_a List of start position and vector of the locus of the first annotation
#
# @param borders_loc_b List of start position and vector of the locus of the second annotation
#
# @see create_vectors()
#
# @return Returns a list describing the matchs/mismatchs for each position of the strings
#
# @remark This function doesn't expect any annotation to be a 'reference'
def old_compare_loci(borders_loc_a, borders_loc_b, debug=False, verbose=False): 
 
    loc_a = create_vectors(borders_loc_a, debug, verbose)
    loc_b = create_vectors(borders_loc_b, debug, verbose)
    
    if debug:
        print(f"loc_a = {loc_a}")
        print(f"loc_b = {loc_b}")
 
    comp_list = [0,0]
    
    # we get the minimum start positions, maximum end positions, and position difference of the loci
    minv = min(loc_a[0], loc_b[0]) # minimum start position
    diff = abs(loc_a[0] - loc_b[0]) # difference between the start positions
    
    if debug:
        print(f"minv = {minv}")
        print(f"diff = {diff}")
    
    # if the two loci don't start at the same position
    if loc_a[0] != loc_b[0]:
        
        if debug:
            print("The two loci do not start at the same position")
        
        # if the locus 'a' starts before the start of locus 'b'
        if minv == loc_a[0]:
            
            if debug:
                print("Locus 'a' start before locus 'b'")
            
            # for every comparison of the numbers at each position in the two strings, we increment by one the match or mismatch value of the return list. We account for the difference in start positions by adding the difference to the locus 'a' codon position retrieval
            for i in range( min(len(loc_a[1]), len(loc_b[1])) + diff ):
                
                # try to know if the current position is in a CDS in the reference or alternative. If we are outside the string of an annotation the value 'False' is assigned
                try:
                    a_in_CDS = loc_a[1][i] in ("1", "2", "3") 
                except IndexError:
                    a_in_CDS = False
                    
                if i >= diff:
                    b_in_CDS = loc_b[1][i-diff] in ("1", "2", "3")
                else:
                    b_in_CDS = False
                
                if debug:
                    print("Range of loop = {min(len(loc_a[1]), len(loc_b[1])) + diff}")
                    print(f"i = {i}")
                    print(f"a_in_CDS = {a_in_CDS}")
                    print(f"b_in_CDS = {b_in_CDS}")
                
                # if we are outside the coordinates of locus 'b' and locus 'a' has a CDS in this position, we increment the mismatch count; else we do nothing
                if i<diff and a_in_CDS:
                    
                    comp_list[0] += 1
                
                # if we are outside the coordinates of locus 'a' and locus 'b' has a CDS in this position, we increment the mismatch count; else we do nothing
                elif i >= len(loc_a[1]) and b_in_CDS:
                    
                    comp_list[0] += 1
                
                # if we are in both coordinates, we increment the match count if both have the same non-zero value at the current position, else we increment the mismatch value if they don't have both '0' as a value
                else:
                    
                    if a_in_CDS and loc_a[1][i] == loc_b[1][i-diff]:
                        
                        comp_list[1] +=1
                
                    elif a_in_CDS or b_in_CDS:
                
                        comp_list[0] += 1
        
        # if the locus'a' starts after the start of locus 'b'
        elif minv == loc_b[0]:
            
            if debug:
                print("Locus 'b' starts before locus 'a'")
            
            # for every comparison of the numbers at each position in the two strings, we increment by one the match or mismatch value of the return list. We account for the difference in start positions by adding the difference to the locus 'b' codon position retrieval
            for i in range( min(len(loc_a[1]), len(loc_b[1])) + diff ):
                
                # try to know if the current position is in a CDS in the reference or alternative. If we are outside the string of an annotation the value 'False' is assigned
                if i >= diff:
                    a_in_CDS = loc_a[1][i-diff] in ("1", "2", "3")
                else:
                    a_in_CDS = False
                    
                try:
                    b_in_CDS = loc_b[1][i] in ("1", "2", "3")
                except IndexError:
                    b_in_CDS = False
                
                if debug:
                    print("Range of loop = {min(len(loc_a[1]), len(loc_b[1])) + diff")
                    print(f"i = {i}")
                    print(f"a_in_CDS = {a_in_CDS}")
                    print(f"b_in_CDS = {b_in_CDS}")
                    
                # if we are outside the coordinates of locus 'a' and locus 'b' has a CDS in this position, we increment the mismatch count; else we do nothing
                if i<diff and b_in_CDS:
                    
                    comp_list[0] += 1
                    
                # if we are outside the coordinates of locus 'b' and locus 'a' has a CDS in this position, we increment the mismatch count; else we do nothing
                elif i >= len(loc_b[1]) and a_in_CDS:
                    
                    comp_list[0] += 1
            
                # if we are in both coordinates, we increment the match count if both have the same non-zero value at the current position, else we increment the mismatch value if they don't have both '0' as a value
                else:
                    
                    if a_in_CDS and loc_a[1][i-diff] == loc_b[1][i]:
                        
                        comp_list[1] +=1
                
                    elif a_in_CDS or b_in_CDS:
                
                        comp_list[0] += 1
    
    # if the two loci start at the same position
    else:
        
        if debug:
            print("The two loci start at the same position")
    
        # for every comparison of the numbers at each position in the two strings, we increment by the match value of the return list of the two string values at this position are equal, or the mismatch value if they are not or if only one annotation has a CDS at this position
        for i in range(len(loc_a[1])):
            
            a_in_CDS = loc_a[1][i] in ("1", "2", "3")
            b_in_CDS = loc_b[1][i] in ("1", "2", "3")
            
            if debug:
                print(f"\nRange of loop = {len(loc_a[1])}")
                print(f"i = {i}")
                print(f"a_in_CDS = {a_in_CDS}")
                print(f"b_in_CDS = {b_in_CDS}")
            
            if a_in_CDS and loc_a[1][i] == loc_b[1][i]:
                          
                comp_list[1] += 1
                
            elif a_in_CDS or b_in_CDS:
                
                comp_list[0] += 1
        
    if verbose:
        print(f"\nResult of the comparison of the locus : {comp_list[1]} matches and {comp_list[0]} mismatches" )
    
    return comp_list

## Main function of this program. Given a reference and alternative path, gets the corresponding GFF files and compares the two annotations to return their structure's identity level
#
# @param ref_path Path of the GFF file describing the reference annotation
#
# @param alt_path Path of the GFF file describing the aternative annotation
#
# @param debug If True, triggers display of many messages intended for debugging the program
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
    
    identities = {}
    
    # for each locus of the reference annotation...
    for locus in ref_annotations:
    
        # if the locus exists in the alternative annotation and if the two loci are on the same strand (same strand number)...
        if locus in alt_annotations and ref_annotations[locus][0] == alt_annotations[locus][0]:
        
            # construct the comparison matrix for the reference and alternative locus annotation...
              
            # using the new version
            if create_strings == False:
                if verbose:
                    print("\nStarting the comparison of the locus " + locus + " using the new program version\n")
                comp_res = compare_loci(ref_annotations[locus][1], alt_annotations[locus][1], debug, verbose)
                
            # using the old version
            else:
                if verbose:
                    print("\nStarting the comparison of the locus " + locus + " using the old program version (structure strings creation)\n")
                comp_res = old_compare_loci(ref_annotations[locus][1], alt_annotations[locus][1], debug, verbose)
            
            # compute identity from the matrix and index it in the identities dictionary
            identities[locus] = round( float(comp_res[1]) / (float(comp_res[1]) + float(comp_res[0]) ) * 100 , 1 )
            
            if verbose:
                print("Finished computing the identity of the annotations of locus " + locus + "\n")
             
        # if the two loci are on different strands or if the locus doesn't exist in the alternative annotation, we assign 0% identity to the locus
        else:
        
            identities[locus] = 0.0
            
            if verbose:
                
                if locus in alt_annotations:
                    print("\nAlternative locus is predicted on a different strand from the reference. Identity : 0%\n")
                else:
                    print("\nLocus is not predicted in alternative annotation. Identity: 0%\n")
                
    return identities


def usage():
    
    # displayed when '-h' or '--help' is given, or when an invalid script call happens
    print("Syntax : path/to/main.py [ -h/--help -d/--debug -v/--verbose -c/--create_strings ] [ -r/--reference <reference_file_path> ] [ -a/--alternative <alternative_file_path> ] ")
    

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
    comparison = annotation_comparison(ref_path, alt_path, debug, verbose, create_strings)
    
    print("\nResult of the comparison of all loci of both annotations :\n")
    pprint.pprint(comparison)
    print("\n\n")
    
if __name__ == "__main__":
    main()
    
#TODO Is the fact that the loci fused by the 'fuse_superloci' function have the same locus ID guaranted ?
#TODO Multiple mRNAs problem (only use the alternative mRNA with the maximum identity ?)
#TODO after running the main program, check if there are loci left in the alternative not present in the reference and assign them 0% identity
#TODO what identity yields comparing overlapping-loci_test with overlapping-loci-alt_test ? (add unit test)
    
    