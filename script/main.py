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
# @return Returns a dictionary of lists containing the number associated with the strand of the locus (1 : direct, 0: reverse) and a list containing the start and end coordinates of all CDS for each locus
def get_gff_borders(path, debug, verbose):
    
    file = open(path, "r") # the file to read
    borders = {} # this variable takes in the borders of each CDS of each gene
    locus_id = "" # the current gene being analyzed
    
    for l in file:
        
        print(l)

        # if we encounter a new gene, we get its ID and create a key in 'borders' with a basic list
        if str(l.split("\t")[2]) == "gene": 
            
            locus_id = l.split("\t")[8].split("")[0].split("\n")[0][3:]
            
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
        
    file.close()
    
    # we return the entire dictionary with all borders
    return borders


## This function retrieves all the CDS coordinates from the given lists of coordinates of the reference (@param ref) and of the alternative (@param alt) and includes them in a unique list of coordinates. The coordinates are sorted in ascending order.
#
# @param ref List of CDS coordinates of the reference annotation returned
#
# @param alt List of CDS coordinates of the alternative annotation returned
#
# @return Returns a list of coordinates compiling all coordinates from both initial lists in ascending order
def get_area_bounds(ref, alt):
    
    bounds=[] # the return list
    i=0
    j=0    
    
    # while we did not reach the end of the coordinates lists...
    while i <= len(ref)-1 and j <= len(alt)-1:
        
        # we add the next coordinate to the result list
        bounds.append(min(ref[i],alt[j]))
        
        if(i==bounds[-1]):
            i+=1
        
        if(j==bounds[-1]):
            j+=1
            
    return bounds

## This function indicates for each couple of bounds in the given list of area bounds (@param area_bounds), if they delimit an area which includes a CDS from the given CDS coordinates list (@param cds_bounds). It is used during the comparison of areas in compare_loci() to know if the reference or alternative have a CDS in the area
#
# @see compare_loci()
#
# @param cds_bounds List of CDS coordinates for an annotation
#
# @param area_bounds List of area bounds
#
# @return Returns a list indicating for each couple of bounds if they include a CDS ('True') or not ('False')
def is_in_cds(cds_bounds, area_bounds):
    
    cds_id=0
    i=0
    in_cds=[] # return list
    
    # while we did not yet reach the end of the bounds list or the CDS coordinates list...
    while (i < len(area_bounds)-1 and cds_id < len(cds_bounds)-1):
        
        current_cds=[cds_bounds[cds_id],cds_bounds[cds_id+1]]
        lb= area_bounds[i]
        ub= area_bounds[i+1]

        # if the current CDS is included in the area bounds, 'True' is used to index it in the return list, else 'False'
        if(current_cds[0]<=lb[0] and lb[1]<=current_cds[1]):
            in_cds.append(True)
        else:
            in_cds.append(False)
            
        # if the current area reaches the end of the current CDS, we skip to the next CDS
        if(ub==current_cds[1]):
            cds_id+=2
            
        i+=1
        
    # if we reached the end of the CDS coordinates list, we continue until all area bounds are exhausted by appending 'False'
    while(i < len(area_bounds)-1):
        in_cds.append(False)
        
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
def compare_loci(ref, alt, debug, verbose):
    
    # we retrieve the bounds of all areas delimited by all the CDS coordinates of both border lists
    area_bounds = get_area_bounds(ref, alt)
    
    # we retrieve the list indicating the areas which include or not a CDS for bot hannotations
    ref_in_CDS = is_in_cds(ref, area_bounds)
    alt_in_CDS = is_in_cds(alt, area_bounds)

    # initialisation of the variables indicating if the annotations have a CDS in the current area, and of the variables keeping track of the current codon position
    ref_in_CDS = True if (ref[0]==area_bounds[0]) else False
    alt_in_CDS = True if (alt[0]==area_bounds[0]) else False
    
    if( ref_in_CDS):
        codon_position_ref=0
    if( alt_in_CDS):
        codon_position_alt=0
  
    
    comparison = [0,0] # return list. First value is mismatches, second value is identities
    bound_id=0
    prev_bounds=area_bounds[0]
    
    # for each comparison area delimited by the bounds in the list 'area_bounds'...
    for bound in area_bounds[1:] :
        
        if(prev_bounds>0):
            
            # if both annotations are in a CDS and have the same codon position, then all codon positions for the rest of the area will be identical, so we add the length of the area to the second value of the comparison list and update both codon positions
            if(ref_in_CDS[bound_id] and alt_in_CDS[bound_id] and codon_position_alt==codon_position_ref):
                comparison[1]+=bound-prev_bounds
                codon_position_alt = (codon_position_alt + (bound-prev_bounds))%3 
                codon_position_ref = (codon_position_ref + (bound-prev_bounds))%3 
                
            # if both annotations are in a CDS but don't have the same codon position, then all codon positions for the rest of the area will be different, so we add the length of the are to the first value of the comparison list and update both codon positions
            elif(ref_in_CDS[bound_id] and alt_in_CDS[bound_id] and codon_position_alt!=codon_position_ref):
                comparison[0]+=bound-prev_bounds
                codon_position_alt = (codon_position_alt + (bound-prev_bounds))%3 
                codon_position_ref = (codon_position_ref + (bound-prev_bounds))%3 
            
            # if only one annotation has a CDS in the comparison area, we add to the first value of the comparison list and update only one codon position
            
            elif(ref_in_CDS[bound_id]):
                comparison[0]+=bound-prev_bounds
                codon_position_ref = (codon_position_ref + (bound-prev_bounds))%3 
                
            
            elif(alt_in_CDS[bound_id]):
                comparison[0]+=bound-prev_bounds
                codon_position_alt = (codon_position_alt + (bound-prev_bounds))%3 
                
            # the case of both annotations being outside of a CDS is not used in the computation of global loci identity, and is ignored
       
        prev_bounds=bound
        
    return comparison


## This function uses the variables created by the function @see compare_loci() to compute the comparison matrix of the comparison of the two areas delimited by lower and upper. It takes a comparison matrix as input to add each area comparison to the global locus comparison matrix.
#
# @see compare_loci()
#
# @param comp_matrix The initial comparison matrix in which to add the result of the area's nucleotide comparisons. It is a list of lists of dimensions (4,4).
#
# @param upper The upper bound on the genome of the locus area to consider
#
# @param lower The lower bound on the genome of the locus area to consider
#
# @param codon_position_ref The codon position of the next CDS nucleotide of the reference annotation
#
# @param codon_position_alt The codon position of the next CDS nucleotide of the alternative annotation
#
# @param ref_in_CDS Boolean indicating if the area includes a reference annotation's CDS
#
# @param alt_in_CDS Boolean indicating if the area includes an alternative annotation's CDS
#
# @param debug If True, triggers display of many messages intended for debugging the program
#
# @param verbose If True, triggers display of more information messages. Default is 'False'
#
# @return Returns the initial comparison matrix (list of lists of dimensions (4,4)) with the area's comparison results added to it
#
def increment_matrix(comp_matrix, upper, lower, codon_position_ref, codon_position_alt, ref_in_CDS, alt_in_CDS, debug, verbose):
    
    # initialisation of the variables indicating in which line (corresponding to the reference) and in which column (corresponding to the alternative) to increment the nucleotide count
    ref_nucl = 0
    alt_nucl = 0
    
    for i in range(0, upper-lower):
        
        # if the area includes a reference annotation's CDS...
        if ref_in_CDS:
            
            # then the line in which to increment the nucleotide count is indicated by the current reference codon position
            ref_nucl = codon_position_ref
            
        else:
            
            # if it does not, then we increment in the line '0'
            ref_nucl = 0
            
        # if the area includes an alternative annotation's CDS...
        if alt_in_CDS:
            
            # then the column in which to increment the nucleotide count is indicated by the current alternative codon position
            alt_nucl = codon_position_alt
            
        else:
            
            # if it does not, then we increment in the column '0'
            alt_nucl = 0
            
        # we increment the corresponding 'cell'
        comp_matrix[ref_nucl][alt_nucl] += 1
        
        # we increment the codon positions if we are in an annotation's CDS
        if ref_in_CDS:
            
            codon_position_ref += 1
            
            if codon_position_ref == 4:
                
                codon_position_ref = 1
                
        if alt_in_CDS:
            
            codon_position_alt += 1
            
            if codon_position_alt == 4:
                
                codon_position_alt = 1
    
    return comp_matrix, codon_position_ref, codon_position_alt


## This function computes the identity of structure strings of two loci from the given comparison 
# matrix (list of lists). Returns the identity as a percentage.
#
# @param matrix The comparison matrix returned by pair_vector_comparison()
#
# @param debug If True, triggers display of many messages intended for debugging the program
#
# @param verbose If True, triggers display of more information messages. Default is 'False'
#
# @see pair_vector_comparison()
#
# @return Returns the identity of the two structure strings describded by the matrix
#
def matrix_to_identity(matrix, debug, verbose):
    
    # number of string positions for which 'caracters' were found to be identical. 'matrix[0][0]' is not taken into account so as to reduce the impact of large introns correctly prédicted
    match = matrix[1][1] + matrix[2][2] + matrix[3][3]
    
    if verbose:
        print("matching codon positions in both annotations = " + str(match))

    # number of string positions for which 'caracters' were found to be identical
    mismatch = matrix[0][1] + matrix[0][2] + matrix[0][3] + matrix[1][0] + matrix[1][2] + matrix[1][3] + matrix[2][0] + matrix[2][1] + matrix[2][3] + matrix[3][0] + matrix[3][1] + matrix[3][2]
    
    if verbose:
        print("mismatching codon positions in both annotations = " + str(mismatch) + "\n")
        
    result = round( float(match) / (float(match) + float(mismatch) ) * 100 , 1 )
    
    if verbose:
        print("Proportion of matching codon positions in both annotations = " + str(result) + "\n")

    return(result)


## Main function of this program. Given a reference and alternative path, gets the corresponding GFF files and compares the two annotations
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
def annotation_comparison(ref_path, alt_path, debug, verbose):

    # get all annotation files and generate the annotation data structure
    ref_annotations = get_gff_borders(ref_path, debug, verbose)
    alt_annotations = get_gff_borders(alt_path, debug, verbose)
    
    identities = {}
    
    # for each locus of the reference annotation...
    for locus in ref_annotations:
    
        # if the two loci are on the same strand (same strand number)...
        if ref_annotations[locus][0] == alt_annotations[locus][0]:
        
            # construct the comparison matrix for the reference and alternative locus annotation
            if verbose:
                print("Creation of comparison matrix of locus " + locus + "\n")
            comp_mat = compare_loci(ref_annotations[locus][1], alt_annotations[locus][1], debug, verbose)
            
            # compute identity from the matrix and index it in the identities dictionary
            ident = matrix_to_identity(comp_mat, debug, verbose)
            identities[locus] = ident
            
            if verbose:
                print("Finished computing the identity of the annotations of locus " + locus + "\n")
             
        # if the two loci are on different strands, we consider they have 0% identity
        else:
        
            identities[locus] = 0.0
            
            if verbose:
                print("Locus predicted on a different strand from the reference in the alternative. Identity : 0%\n")
            
    return identities


def usage():
    
    # displayed when '-h' or '--help' is given, or when an invalid script call happens
    print("Syntax : path/to/main.py [ -h/--help -d/--debug -v/--verbose ] [ -r/--reference <reference_file_path> ] [ -a/--alternative <alternative_file_path> ] ")
    

def main():
    
    # we retrieve all script call options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hdvcr:a:", ["help", "debug", "verbose", "chatty", "reference=", "alternative="])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)
        
    # initialisation of the display parameters
    debug = False
    verbose = False
    
    # we get the values given for each parameter
    for o, a in opts:
        if o in ("-d", "--debug"):
            debug = True
        elif o in ("-v", "--verbose"):
            verbose = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-r", "--reference"):
            ref_path = a
        elif o in ("-a", "--alternative"):
            alt_path = a
        else:
            assert False, "unhandled option"
            
    # call of the annotation_comparison function
    comparison = annotation_comparison(ref_path, alt_path, debug, verbose)
    
    pprint.pprint(comparison)
    
if __name__ == "__main__":
    main()
    
    
    
    