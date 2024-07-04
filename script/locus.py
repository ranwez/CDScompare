#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 21:04:04 2024

@author: vetea
"""


## This class represents an annotation's locus identified from reading the
# associated GFF file with the get_gff_borders function. It posesses multiple 
# attributes describing the locus
#
# @see get_gff_borders()
class Locus:
    
    ## This method initialises the class with default values for each attribute
    #
    # @param name (optional) ID of the locus
    #
    # @param mRNAs (optional) Dictionary listing all mRNAs of the locus. Each 
    # key corresponds to the ID of the mRNA and each value to the list of all 
    # CDS coordinates as retrieved by the get_gff_borders function
    #
    # @param start (optional) Start coordinate of the locus
    #
    # @param end (optional) End coordinate of the locus
    #
    # @param direction (optional) String indicating if the locus is on the 
    # 'direct' or 'reverse' strand
    #
    # @see get_gff_borders()
    def __init__(self, name="", mRNAs=None, start=-1, end=-1, direction=""):
        if not mRNAs:
            mRNAs = {}
        
        ### name (identifier) of the locus
        self.name = name
        
        ### Dictionary of lists containing the CDS list of each mRNA of the locus
        self.mRNAs = mRNAs.copy()
        
        ### start coordinate of the locus on the sequence
        self.start = start
        
        ### end coordinate of the locus on the sequence
        self.end = end
        
        ### string indicating if the locus is on direct or reverse strand
        self.direction = direction
        
    ## This method is used to retrieve the 'mRNAs' attribute of an instance
    #
    # @returns Returns the 'mRNAs' attribute of the class instance
    def mRNAs(self):
        return self.mRNAs
    
    ## Sets the value of the 'mRNAs' attribute of the class instance
    #
    # @param value Value to be assigned to the 'mRNAs' attribute
    def set_mRNAs(self, value):
        self.mRNAs = value
    
    ## This method is used to verify if all mRNAs of a dictionary are present 
    # in the class instance's 'mRNAs' attribute
    #
    # @param **mrnas Dictionary with mRNA names as keys and CDS coordinates
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
        
        
    ## Reverses the coordinates of all mRNAs of the class instance
    #
    # @param cluster_end End position coordinate of the cluster containing
    # the locus
    #
    # @remark This method doesn't return anything and modifies directly
    # the class instance. This method is used in case of loci on the reverse
    # strand, and transforms the coordinates into coordinates from the end
    # of the parent cluster
    def reverse(self, cluster_end):
        new_mRNAs = {}
        for mRNA_id, mRNA in self.mRNAs.items():
            new_mRNAs[mRNA_id] = []
            for i in range(len(mRNA)-1, -1, -1):
                new_mRNAs[mRNA_id].append(cluster_end-mRNA[i])
        self.set_mRNAs(new_mRNAs)


    ## This function was added for convenience as a way to easily retrieve the
    # attribute values of the class instance as a formatted string
    #
    # @returns Returns a string describing the values of each attribute of the
    # class instance
    def show_init(self):
        return f"Locus(name='{self.name}', mRNAs={self.mRNAs}, start={self.start}, end={self.end}, direction='{self.direction}'"