#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 11:57:29 2024

@author: vetea
"""

debug = True

path = "../data/tests/test_basique.gff3"

# This function expects a string corresponding to the file path of the GFF file to read, and returns
# a dictionary of lists with the keys corresponding to a locus identifier, and the values corresponding to 
# a list of the start and end position of each coding sequence ('CDS') of the locus 
def get_gff_borders(path):
    
    file = open(path, "r")
    borders = {} # This variable takes in the borders of each CDS of each gene
    locus_id = "" # the current gene being analyzed
    
    for l in file:

        # if we encounter a new gene, we get its ID and create a key in 'borders' with an empty list
        if str(l.split("\t")[2]) == "gene": 
            
            locus_id = l.split("\t")[8].split(";")[0][3:]
            
            borders[locus_id] = []     
            
            if debug :
                print("Reading the locus " + locus_id)

        # if we encounter a CDS, we add its start and end positions to corresponding gene key in 'borders'
        if str(l.split("\t")[2]) == "CDS":
            
            borders[locus_id].append(l.split("\t")[3])
            borders[locus_id].append(l.split("\t")[4])
            
            if debug:
                print("Adding borders to " + locus_id + " : " + l.split("\t")[3] + ", " + l.split("\t")[4])
    
    file.close()
    
    # we return the entire dictionary with all borders
    return borders

