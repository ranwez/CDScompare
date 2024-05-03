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
                print("Reading the locus " + locus_id)

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
# same keys, and as values a string describing the codon position of each nucleotide (1,2,3, or 
# 0 in the case of a non-CDS nucleotide) of the locus/gene
#
# @param borders The dictionary containing all start-end coordinates of the annotation's CDS
#
# @see get_gff_borders()
#
# @return Returns a dictionary of strings describing the annotation structure of each locus
def create_vectors(borders):

    verbose = True # temporary workaround so that unitary tests work
    
    vectors = {} # this variable takes in the strings of gene annotation structure for each gene
    
    # for each locus indexed in the 'borders' dictionary, we create a new key in 'vectors' and 
    # create the structure string as its value
    for gene in borders:
        
        if verbose:
            print("\nConverting the locus " + gene + " with coordinates " + str(borders[gene]))
        
        vectors[gene] = ""
        
        # this variable takes the codon position of the next CDS nucleotide and loops
        # between the values 1, 2, and 3
        codon_pos = 1 
        
        # this variable indicates if we are in an exon/CDS or not.
        # it is incremented at each transition between annotations (each element in the 'borders' list)
        # to represent the exon-intron change along the gene
        in_exon = 1
        
        # for each coordinate indexed for this locus in 'borders'
        for i in range(len(borders[gene])-1):
        
            # if we are in a CDS, we append the numbers 1, 2, and 3 (with looping) 
            if in_exon % 2 == 1:
                
                # for each nucleotide between this coordinate and the next...
                for j in range( borders[gene][i+1] - borders[gene][i] ):
                    
                    # we append the codon position to the structure string 
                    vectors[gene] += str(codon_pos)
                    
                    # we increment the codon position with looping
                    if codon_pos == 3:
                        codon_pos = 1
                    else:
                        codon_pos += 1
                
                
            # if we are not in a CDS, we append 0
            if in_exon % 2 == 0:
                
                # for each nucleotide between this coordinate and the next...
                for j in range( borders[gene][i+1] - borders[gene][i] ):                
                
                    vectors[gene] += "0"
            
            in_exon += 1    
            
        if verbose:
            print("\nStructure string of the " + gene + " locus :\n" + vectors[gene] + "\n")
            
    return vectors
  
## This function expects two structure strings corresponding to two annotations of the same
# genome, and returns a list of lists describing the comparison of each position of each string
#
# @param vect_a Structure string of the first annotation
#
# @param vect_b Structure string of the second annotation
#
# @see create_vectors()
#
# @return Returns a list of lists (matrix) describing the matchs/mismatchs for each position
# of the strings
#
# @remark This function expects both annotations to be of the same size, but doesn't expect any to be a # 'reference'
def pair_vector_comparison(vect_a, vect_b):
 
    verbose = True # temporary workaround so that unitary tests work    
 
    # initialising the return matrix with 4 'lines'  and 4 'columns' (for 0, 1, 2, 3)
    comp_matrix = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    
    # for every comparison of the numbers at each position in the two strings, we increment by one 
    # the corresponding 'cell'
    for i in range(len(vect_a)):
        
        comp_matrix[int(vect_a[i])][int(vect_b[i])] += 1
        
    if verbose:
        print("\n" + str(comp_matrix[0]) + "\n" + str(comp_matrix[1]) + "\n" + str(comp_matrix[2]) + "\n" + str(comp_matrix[3]) + "\n" )
    
    return comp_matrix


def usage():
    
    print("Syntax : path/to/main.py [ -h/--help -d/--debug -v/--verbose ] [ -i/--input <input_folder_path> ]")
    

def main():
    
    # get all script call options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hdvi:", ["help", "debug", "verbose", "input="])
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
            path = a
        else:
            assert False, "unhandled option"
            
    path2 = "../data/tests/identical_test.gff3"
    path3 = "../data/tests/minus-CDS_test.gff3"
    path4 = "../data/tests/fusion_test.gff3"
    path5 = "../data/tests/shift_test.gff3"
            
    vect = create_vectors( get_gff_borders(path) )
    vect2 = create_vectors( get_gff_borders(path2) )
    vect3 = create_vectors( get_gff_borders(path3) )
    vect4 = create_vectors( get_gff_borders(path4) )
    vect5 = create_vectors( get_gff_borders(path5) )
    
    pair_vector_comparison(vect["chr2A_00611930"], vect2["chr2A_00611930"])
    pair_vector_comparison(vect["chr2A_00611930"], vect3["chr2A_00611930"])
    pair_vector_comparison(vect["chr2A_00611930"], vect4["chr2A_00611930"])
    pair_vector_comparison(vect["chr2A_00611930"], vect5["chr2A_00611930"])
    
if __name__ == "__main__":
    main()
    
    
    
    