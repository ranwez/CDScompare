#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 21:04:24 2024

@author: vetea
"""


#TODO update docu
## This class represents clusters of overlapping loci computed by the 
# function construct_clusters. It posesses only one attribute, 'clusters',
# which is a dictionary of dictionaries of instances of the class 'Locus'
# indicating the loci attributed to each cluster.
#
# @see construct_clusters()
#
# @see Locus
class Cluster:
    
    ## This method initialises the class with an empty dictionary
    def __init__(self, name, loci=None):
        self.name = name
        if not loci:
            loci = {"ref": [], "alt": []}
        self.loci = loci.copy()
        
    ## This method is used to retrieve the 'clusters' attribute of an instance
    #
    # @returns Returns the 'clusters' attribute of the class instance
    def get_loci(self):
        return self.loci.copy()
    
    #TODO docu
    def append_to_loci(self, annotation, value):
        self.loci[annotation].append(value)
        
    #TODO docu
    def get_end(self):
        return self.end
    
    #TODO docu
    def set_end(self, value):
        self.end = value
    
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
    
    def get_details(self):
        result = {'ref': [], 'alt': []}
        for loc in self.get_loci()['ref']:
            result['ref'].append(loc.mRNAs)
        for loc in self.get_loci()['alt']:
            result['alt'].append(loc.mRNAs)
        return result
        
        
        
        