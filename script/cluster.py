#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 21:04:24 2024

@author: vetea
"""


## This class represents clusters of overlapping loci computed by the 
# function construct_clusters. It posesses only three attributes: 'loci',
# which is a dictionary of lists of instances of the class 'Locus'
# indicating the loci attributed to each cluster; 'name'; and 'end' which is 
# the end coordinate of the last loci of the cluster
#
# @see construct_clusters()
#
# @see Locus
class Cluster:
    
    ## This method initialises the class
    #
    # @param name (optional) the name of the cluster
    #
    # @param loci (optional) the lists of loci for the cluster of the reference
    # and alternative annotations, as a dictionary of lists
    #
    # @param end (optional) the end coordinate of the last locus of the cluster
    def __init__(self, name, loci=None, end=-1):
        
        ### name of the cluster
        self.name = name
        if not loci:
            loci = {"ref": [], "alt": []}
        
        ### list of Locus class instances present in the cluster
        self.loci = loci.copy()
        
        ### end coordinate of the cluster. Equal to the end coordinate of 
        # the last cluster's locus
        self.end = end
        
    ## Retrieves the 'loci' attribute of the instance
    #
    # @returns Returns the 'loci' attribute of the class instance
    def get_loci(self):
        return self.loci.copy()
    
    ## Appends the given value to the list of the given annotation in the 
    # 'loci' attribute
    #
    # @param annotation string ('ref' or 'alt') indicating the origin of the
    # value
    #
    # @param value Value to append to the 'loci' attribute
    def append_to_loci(self, annotation, value):
        self.loci[annotation].append(value)
        
    ## Retrieves the 'end' attribute of the instance
    #
    # @returns Returns the 'end' attribute of the class instance
    def get_end(self):
        return self.end
    
    ## Sets the given value as the 'end' attribute
    #
    # @param value Value to give to the 'end' attribute of the class instance
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
    
    ## Retrieves the class instance's 'loci' attribute as a human-readable 
    # dictionary of lists detailing each locus' mRNAs
    #
    # @returns Returns a dictionary of list detailing the content of the 
    # instance's 'loci' attribute
    def get_details(self):
        result = {'ref': [], 'alt': []}
        for loc in self.get_loci()['ref']:
            result['ref'].append(loc.mRNAs)
        for loc in self.get_loci()['alt']:
            result['alt'].append(loc.mRNAs)
        return result
        
        
        
        