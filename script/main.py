#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 11:57:29 2024

@author: vetea
"""

##@package main
# This script is used to compute the distance between multiple structural annotations of a same genome. # It expects as input the path to the folder containing all annotation files (in GFF format), 
# displays the computed distances between all annotation pairs, and returns a dictionary of lists of
# lists detailing the matchs/mismatchs between the two annotations' structure string

import getopt
import sys
import os
import pprint


## This function expects a string corresponding to the file path of the GFF file to read, and returns
# a dictionary of lists with the keys corresponding to a locus identifier, and the values 
# corresponding to a list of the start and end position of each coding sequence ('CDS') of the locus 
#
# @param path Path of the file to read
#
# @return Returns a dictionary of lists containing the start and end coordinates of 
# all CDS of each locus
def get_gff_borders(path):
    
    verbose = True # temporary workaround so that unitary tests work
    
    file = open(path, "r") # the file to read
    borders = {} # this variable takes in the borders of each CDS of each gene
    locus_id = "" # the current gene being analyzed
    
    for l in file:

        # if we encounter a new gene, we get its ID and create a key in 'borders' with an empty list
        if str(l.split("\t")[2]) == "gene": 
            
            locus_id = l.split("\t")[8].split(";")[0].split("\n")[0][3:]
            
            borders[locus_id] = []     
            
            if verbose :
                print("\nReading the locus " + locus_id)

        # if we encounter a CDS, we add its start and end positions to corresponding gene key in 'borders'
        if str(l.split("\t")[2]) == "CDS":
            
            borders[locus_id].append(int(l.split("\t")[3]))
            borders[locus_id].append(int(l.split("\t")[4]))
            
            if verbose:
                print("Adding borders to " + locus_id + " : " + l.split("\t")[3] + ", " + l.split("\t")[4])
    
    file.close()
    
    # we return the entire dictionary with all borders
    return borders


## This function expects a dictionary with keys corresponding to gene IDs and values corresponding
# to a list of all CDS coordinates (start and end) of the gene. It returns a dictionary with the
# same keys, and as values an int indicating the start of the locus and a string describing the codon
# position of each nucleotide (1,2,3, or 0 in the case of a non-CDS nucleotide) of the locus/gene
#
# @param borders The dictionary containing all start-end coordinates of the annotation's CDS
#
# @see get_gff_borders()
#
# @return Returns a dictionary of lists describing the start of the locus and the annotation structure of each locus
def create_vectors(borders):

    verbose = True # temporary workaround so that unitary tests work
    
    vectors = {} # this variable takes in the strings of gene annotation structure for each gene
    
    # for each locus indexed in the 'borders' dictionary, we create a new key in 'vectors' and 
    # create the structure string as its value
    for gene in borders:
        
        if verbose:
            print("\nConverting the locus " + gene + " with coordinates " + str(borders[gene]))
        
        # if the locus is on the direct strand
        if borders[gene][1] > borders[gene][0]:
            
            # we get the start of the locus
            vectors[gene] = [ borders[gene][0] , "" ]
            
        # if the locus is on the reverse strand
        else:
            # we get the start (the last CDS value since the locus is reversed) of the locus
            vectors[gene] = [ borders[gene][-1] , "" ]
        
        # this variable takes the codon position of the next CDS nucleotide and loops
        # between the values 1, 2, and 3
        codon_pos = 1 
        
        # this variable indicates if we are in an exon/CDS or not.
        # it is incremented at each transition between annotations (each element in the 'borders' list)
        # to represent the exon-intron change along the gene
        in_exon = 1
        
        if borders[gene][1] > borders[gene][0]:
        
            # for each coordinate indexed for this locus in 'borders'
            for i in range(len(borders[gene])-1):
            
                # if we are in a CDS, we append the numbers 1, 2, and 3 (with looping) 
                if in_exon % 2 == 1:
                    
                    # for each nucleotide between this coordinate and the next...
                    for j in range( borders[gene][i+1] - borders[gene][i] ):
                        
                        # we append the codon position to the structure string 
                        vectors[gene][1] += str(codon_pos)
                        
                        # we increment the codon position with looping
                        if codon_pos == 3:
                            codon_pos = 1
                        else:
                            codon_pos += 1
                    
                    
                # if we are not in a CDS, we append 0
                if in_exon % 2 == 0:
                    
                    # for each nucleotide between this coordinate and the next...
                    for j in range( borders[gene][i+1] - borders[gene][i] ):                
                    
                        vectors[gene][1] += "0"
                
                in_exon += 1    
            
        elif borders[gene][1] < borders[gene][0]:

            # for each coordinate indexed for this locus in 'borders' in reverse order
            for i in range(len(borders[gene])-1, 0, -1):
            
                # if we are in a CDS, we append the numbers 1, 2, and 3 (with looping) 
                if in_exon % 2 == 1:
                    
                    # for each nucleotide between this coordinate and the next...
                    for j in range( borders[gene][i-1] - borders[gene][i] ):
                        
                        # we append the codon position to the structure string 
                        vectors[gene][1] += str(codon_pos)
                        
                        # we increment the codon position with looping
                        if codon_pos == 3:
                            codon_pos = 1
                        else:
                            codon_pos += 1
                    
                    
                # if we are not in a CDS, we append 0
                if in_exon % 2 == 0:
                    
                    # for each nucleotide between this coordinate and the next...
                    for j in range( borders[gene][i-1] - borders[gene][i] ):
                    
                        vectors[gene][1] += "0"
                
                in_exon += 1   
        
        if verbose:
            print("\nStructure string of the " + gene + " locus :\n" + vectors[gene][1] + "\n")
            
    return vectors


## This function creates a dictionary of dictionaries of structure strings corresponding to the
# structure string for each locus of each annotation found in the given list of file paths.
#
# @param file_list List of file paths corresponding to all annotations analyzed (including reference)
#
# @return Returns a dictionary of dictionaries of the structure strings for each locus of 
# each annotation
#
def create_all_vectors(file_list):
    
    verbose = True # temporary workaround so that unitary tests work    
    
    annotations = {}
    
    for file in file_list:
        
        filename = file.split("/")[-1]
        
        annotations[filename] = ( create_vectors( get_gff_borders(file) ) )
        
    if verbose:
        print("Finished creating all annotation structure strings")
        
    return annotations
  
    
## This function expects two structure strings corresponding to two annotations of the same
# genome, and returns a list of lists describing the comparison of each position of each string
#
# @param loc_a List of start position and vector of the locus of the first annotation
#
# @param loc_b List of start position and vector of the locus of the second annotation
#
# @see create_vectors()
#
# @return Returns a list of lists (matrix) describing the matchs/mismatchs for each position
# of the strings
#
# @remark This function expects both annotations to be of the same size, but doesn't expect any to be a # 'reference'
def pair_vector_comparison(loc_a, loc_b):
 
    verbose = True # temporary workaround so that unitary tests work    
 
    # initialising the return matrix with 4 'lines'  and 4 'columns' (for 0, 1, 2, 3)
    comp_matrix = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    
    # we get the minimum start positions, maximum end positions, and position difference of the loci
    minv = min(loc_a[0], loc_b[0]) # minimum start position
    diff = abs(loc_a[0] - loc_b[0]) # difference between the start positions
    
    # if the two loci don't start at the same position
    if loc_a[0] != loc_b[0]:
        
        # if the locus 'a' starts before the start of locus 'b'
        if minv == loc_a[0]:
            
            # for every comparison of the numbers at each position in the two strings, we increment by one the corresponding 'cell'. We account for the difference in start positions by adding the difference to the locus 'a' codon position retrieval
            for i in range( min(len(loc_a[1]), len(loc_b[1])) + diff ):
                
                # if we are outside the coordinates of locus 'b', we replace its codon position with 0
                if i<diff:
                    
                    comp_matrix[ int(loc_a[1][i]) ][0] += 1
                
                # if we are outside the coordinates of locus 'a', we replace its codon position with 0
                elif i >= len(loc_a[1]):
                    
                    comp_matrix[0][ int(loc_b[1][i-diff]) ] += 1
                
                # if we are in both coordinates, get each codon position
                else:
                
                    comp_matrix[ int(loc_a[1][i]) ][ int(loc_b[1][i-diff]) ] += 1
        
        # if the locus'a' starts after the start of locus 'b'
        elif minv == loc_b[0]:
            
            # for every comparison of the numbers at each position in the two strings, we increment by one the corresponding 'cell'. We account for the difference in start positions by adding the difference to the locus 'b' codon position retrieval
            for i in range( min(len(loc_a[1]), len(loc_b[1])) + diff ):
                
                # if we are outside the coordinates of locus 'a', we replace its codon position with 0
                if i<diff:
                    
                    comp_matrix[0][ int(loc_b[1][i]) ] += 1
                    
                # if we are outside the coordinates of locus 'b', we replace its codon position with 0
                elif i >= len(loc_b[1]):
                    
                    comp_matrix[ int(loc_a[1][i-diff]) ][0] += 1
            
                # if we are in both coordinates, get each codon position
                else:
                    
                    comp_matrix[ int(loc_a[1][i-diff]) ][ int(loc_b[1][i]) ] += 1
    
    # if the two loci start at the same position
    else:
    
        # for every comparison of the numbers at each position in the two strings, we increment by one 
        # the corresponding 'cell'
        for i in range(len(loc_a[1])):
            
            comp_matrix[ int(loc_a[1][i]) ][ int(loc_b[1][i]) ] += 1
        
    if verbose:
        print("\n" + str(comp_matrix[0]) + "\n" + str(comp_matrix[1]) + "\n" + str(comp_matrix[2]) + "\n" + str(comp_matrix[3]) + "\n" )
    
    return comp_matrix


## This function computes the identity of structure strings of two loci from the given comparison 
# matrix (list of lists). Returns the identity as a percentage.
#
# @param matrix The comparison matrix returned by pair_vector_comparison()
#
# @see pair_vector_comparison()
#
# @return Returns the identity of the two structure strings describded by the matrix
#
def matrix_to_identity(matrix):
    
    verbose = True # temporary workaround so that unitary tests work
    
    # number of string positions for which 'caracters' were found to be identical. 'matrix[0][0]' is not taken into account so as to reduce the impact of large introns correctly prédicted
    match = matrix[1][1] + matrix[2][2] + matrix[3][3]
    
    if verbose:
        print("matching codon positions in both annotations = " + str(match))

    # number of string positions for which 'caracters' were found to be identical
    mismatch = matrix[0][1] + matrix[0][2] + matrix[0][3] + matrix[1][0] + matrix[1][2] + matrix[1][3] + matrix[2][0] + matrix[2][1] + matrix[2][3] + matrix[3][0] + matrix[3][1] + matrix[3][2]
    
    if verbose:
        print("mismatching codon positions in both annotations = " + str(mismatch) + "\n")

    return( round( float(match) / (float(match) + float(mismatch) ) * 100 , 1 ) )


## This function gets all GFF files from the given folder path and returns a list of file paths
# corresponding to each GFF file encountered
#
# @param folder_path Path of the folder from which to retrieve GFF files
#
def get_files(folder_path):
    
    verbose = True # temporary workaround so that unitary tests work
    
    file_list = []
    
    # code adapted from https://stackoverflow.com/questions/10377998/how-can-i-iterate-over-files-in-a-given-directory
    folder = os.fsencode(folder_path)

    # for each file in the given folder
    for file in os.listdir(folder):
        
        filename = os.fsdecode(file)
        
        # if the file has a GFF/GTF extension, we add its path to the list
        if filename.endswith(".gff") or filename.endswith(".gff3") or filename.endswith(".gtf"):
            file_list.append(folder_path+filename)
            continue
        else:
            continue
        
    if verbose:
        
        print(str(file_list))

    return file_list


## Main function of this program. Given a folder path, gets all GFF files in it, 
# creates structure strings dictionaries for each annotation, and compares each annotation with
# the indicated reference annotation.
#
# @param folder_path Path of the folder from which to get the GFF files
#
# @param ref_name Complete name of the reference annotation file (without the folder path)
#
# @return Returns a dictionary of dictionaries of floats corresponding to the structure
# string identity between each locus of each annotation compared to those of the reference
#
# @remark Loci found in one annotation but not the other are ignored
def annotation_comparison(folder_path, ref_name):
    
    verbose = True # temporary workaround so that unitary tests work    

    # get all annotation files and generate the annotation data structure
    files = get_files(folder_path)
    annotations = create_all_vectors(files)
    
    identities = {}

    # for each annotation...
    for ann in annotations:

        identities[ann] = {}
        
        # for each locus of the annotation...
        for locus in annotations[ann]:
            
            # if the locus is present in the reference annotation...
            if locus in annotations[ref_name]:
                
                # construct the comparison matrix
                if verbose:
                    print("Comparison matrix of locus " + locus + " of annotation " + str(ann) + "with equivalent of reference annotation")
                comp_mat = pair_vector_comparison( annotations[ref_name][locus], annotations[ann][locus])
                
                # compute identity from the matrix and index it in the identities dictionary
                ident = matrix_to_identity(comp_mat)
                identities[ann][locus] = ident
                
                if verbose:
                    print("identity of comparison of locus " + locus + " of annotation " + str(ann) + "with equivalent of reference annotation : " + str(ident) + "\n")
                
            
    return identities


def compare(ref, alt):
    
    i = 0
    j = 0
    cpa = 1
    cpb = 1
    upper = 0
    lower = 0
    
    while i <= len(ref)-1 and j <= len(alt)-1:
        
        if i == len(ref)-1:
            
            for k in range(j, len(alt)-1):
                
                lower = alt[k]
                upper = alt[k+1]
                
                print(lower, upper)
                
            i += 1
            
        elif j == len(alt)-1:
            
            for k in range(i, len(ref)-1):
                
                lower = ref[k]
                upper = ref[k+1]
                
                print(lower, upper)
                
            j += 1

        elif ref[i] == alt[j]:
    
            if ref[i+1] <= alt[j+1]:
                
                i += 1
                
            else:
                
                j += 1
    
        else:
            
            if min(ref[i], alt[j]) == ref[i]:
                
                lower = ref[i]
                i += 1
                
            else:
                
                lower = alt[j]
                j += 1
            
            if i == len(ref) or j == len(alt): 
            
                upper = max(ref[-1], alt[-1]) 
            
            else:
                
                upper = min(ref[i], alt[j])
            
            print(lower, upper)

def usage():
    
    print("Syntax : path/to/main.py [ -h/--help -d/--debug -v/--verbose ] [ -i/--input <input_folder_path> ] [ -r/--reference <reference_file_name> ]")
    

def main():
    
    # get all script call options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hdvi:r:", ["help", "debug", "verbose", "input=", "reference="])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)
        
    global debug
    global verbose
    global test
    debug = False
    verbose = False
    
    for o, a in opts:
        if o in ("-d", "--debug"):
            debug = True
        elif o in ("-v", "--verbose"):
            verbose = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-i", "--input"):
            folder_path = a
        elif o in ("-r", "--reference"):
            ref_name = a
        else:
            assert False, "unhandled option"
            
    print(ref_name)
            
    comparison = annotation_comparison(folder_path, ref_name)
    
    pprint.pprint(comparison)
    
if __name__ == "__main__":
    main()
    
    
    
    