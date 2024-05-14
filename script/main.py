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
            
            locus_id = l.split("\t")[8].split(";")[0].split("\n")[0][3:]
            
            # the locus list is intialised with a strand number of '1', which is the corrected if necessary
            borders[locus_id] = [1, []]
            
            if verbose :
                print("\nReading the locus " + locus_id)

        # if we encounter a CDS, we add its start and end positions to corresponding gene key in 'borders'
        if str(l.split("\t")[2]) == "CDS":
            
            borders[locus_id[1]].append(int(l.split("\t")[3]))
            borders[locus_id[1]].append(int(l.split("\t")[4]))
            
            if verbose:
                print("Adding borders to " + locus_id + " : " + l.split("\t")[3] + ", " + l.split("\t")[4])
            # if the retrieved borders indicate the locus is on the reverse strand, we change the strand number to '0'
            if borders[locus_id[0]] == 1 and borders[locus_id[1][1]] < borders[locus_id[1][0]]:
                
                borders[locus_id[0]] = 0
        
    file.close()
    
    # we return the entire dictionary with all borders
    return borders


## This function compares two annotations' loci from their CDS border list returned by the function get_gff_borders and creates a comparison matrix detailing the identities and differences between the two annotations's codon position structure
#
# @see get_gff_borders
#
# @param ref The reference annotation's border list
#
# @param alt The alternative annotation's border list
#
# @param debug If True, triggers display of many messages intended for debugging the program
#
# @param verbose If True, triggers display of more information messages. Default is 'False'
#
# @return Returns a list of list (matrix) indicating what codon position is indicated in the reference and alternative annotations (values : 1, 2, 3 (CDS), or 0 (intron))
#
def compare_loci(ref, alt, debug, verbose):
    
    # initialisation of the border list iterators (i iterator of the reference annotation and j iterator of the alternative annotation)
    i = 0
    j = 0
    
    # initialisation of the variables indicating the codon position of the next CDS nucleotide
    codon_position_ref = 1
    codon_position_alt = 1
    
    # initialisation of the comparison area delimitation variables. 'upper' indicates the upper border of the zone of comparison of the two annotations and 'lower' indicates the lower border
    upper = 0
    lower = 0
    
    # initialisation of the comparison matrix to return
    comp_matrix = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    
        
    ref_in_CDS = True if (ref[0]==min(ref[0],alt[0])) else False
    alt_in_CDS = True if (alt[0]==min(ref[0],alt[0])) else False

   
    
    # while the comparison is not outside of the bounds of the border lists...
    while i <= len(ref)-1 and j <= len(alt)-1:
        
        # if the comparison reached the last border coordinate of the reference list...
        if i == len(ref)-1:
            
            # we add the nucleotide comparisons of the last comparison area of the reference annotation (which is inside a CDS)
            
            ref_in_CDS = True
            
            lower = alt[j]
            upper = alt[j+1]
            alt_in_CDS = not alt_in_CDS
                
            # we call the increment_matrix function to add the area's nucleotide comparisons to the locus' comparison matrix
            comp_matrix, codon_position_ref, codon_position_alt = increment_matrix(comp_matrix, upper, lower, codon_position_ref, codon_position_alt, ref_in_CDS, alt_in_CDS, debug, verbose)
            
            # then we add the nucleotide comparisons for all the comparison areas after the end of the reference annotation (outside a reference CDS)
            
            ref_in_CDS = False
            
            for k in range(j+1, len(alt)-1):
                
                lower = alt[k]
                upper = alt[k+1]
                alt_in_CDS = not alt_in_CDS
                
                # we call the increment_matrix function to add the area's nucleotide comparisons to the locus' comparison matrix
                comp_matrix, codon_position_ref, codon_position_alt = increment_matrix(comp_matrix, upper, lower, codon_position_ref, codon_position_alt, ref_in_CDS, alt_in_CDS, debug, verbose)
                
            i += 1
            
        # if the comparison reached the last border coordinate of the alternative list...            
        elif j == len(alt)-1:
            
            # we add the nucleotide comparisons of the last comparison area of the alternative annotation (which is inside a CDS)            
            
            alt_in_CDS = True
                
            lower = ref[i]
            upper = ref[i+1]
            ref_in_CDS = not ref_in_CDS
                
            # we call the increment_matrix function to add the area's nucleotide comparisons to the locus' comparison matrix
            comp_matrix, codon_position_ref, codon_position_alt = increment_matrix(comp_matrix, upper, lower, codon_position_ref, codon_position_alt, ref_in_CDS, alt_in_CDS, debug, verbose)
            
            # then we add the nucleotide comparisons for all the comparison areas after the end of the alternative annotation (outside a reference CDS)
            
            alt_in_CDS = False
            
            for k in range(i+1, len(ref)-1):
                
                lower = ref[k]
                upper = ref[k+1]
                ref_in_CDS = not ref_in_CDS
                
                # we call the increment_matrix function to add the area's nucleotide comparisons to the locus' comparison matrix
                comp_matrix, codon_position_ref, codon_position_alt = increment_matrix(comp_matrix, upper, lower, codon_position_ref, codon_position_alt, ref_in_CDS, alt_in_CDS, debug, verbose)
                
            j += 1

        # if the genomic positions indicated by i and j are the same for both annotations...
        elif ref[i] == alt[j]:
    
            # then we increment i or j to circumvent a problem with the function counting two times the same area
            
            if ref[i+1] <= alt[j+1]:
                
                i += 1
                
            else:
                
                j += 1
    
        # if the coordinates i and j do not indicate the same genomic position and we did not reach the end of the border lists...
        else:
            
            # if the position indicated by i is less than the position of j...
            if min(ref[i], alt[j]) == ref[i]:
                
                # then the position of i is used as the lower border...
                lower = ref[i]
                
                # and we can deduce that the comparison area will include a reference annotation's CDS
                ref_in_CDS = True if i%2 == 0 else False
                
                i += 1
                
            # if the position indicated by j is less than the position of i...
            else:
                
                # then the position of j is used as the lower border...
                lower = alt[j]
                
                # and we can deduce that the comparison area will include an alternative annotation's CDS
                alt_in_CDS = True if j%2 == 0 else False
                
                j += 1
            
            # if we are at a last position of the border lists...
            if i == len(ref) or j == len(alt): 
            
                # if the reference's last position is greater than the alternative's...
                if ref[-1] > alt[-1]:
                    
                    # then we use it as an upper border of the comparison area
                    upper = ref[-1]
                    
                    # and we can deduce that the comparison area will include a reference annotation's CDS
                    ref_in_CDS = True

                # if the alternative's last position is greater than the reference's...                    
                else:
                    
                    # then we use it as an upper border of the comparison area
                    upper = alt[-1]
                    
                    # and we can deduce that the comparison area will include an alternative annotation's CDS
                    alt_in_CDS = True
            
            # if we are not at a last position (still inside both loci)...
            else:
                
                # if the next position in both lists is in the reference...
                if ref[i] < alt[j]:
                    
                    # then we use the next reference annotation position as an upper comparison area border...
                    upper = ref[i]
                    
                    # and we can deduce that the comparison area will include a reference annotation's CDS
                    ref_in_CDS = True if i%2 == 1 else False
                    
                # if the next position in both lists is in the reference...
                else:
                    
                    # then we use the next alternative annotation position as an upper comparison area border...
                    upper = alt[j]
                    
                    # and we can deduce that the comaprison area will include an alternative annotation's CDS
                    alt_in_CDS = True if j%2 == 1 else False
                    
            # we call the increment_matrix function to add the area's nucleotide comparisons to the locus' comparison matrix
            comp_matrix, codon_position_ref, codon_position_alt = increment_matrix(comp_matrix, upper, lower, codon_position_ref, codon_position_alt, ref_in_CDS, alt_in_CDS, debug, verbose)
            
    return comp_matrix


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
    
    
    
    