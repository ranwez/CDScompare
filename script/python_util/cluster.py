#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 21:04:24 2024

@author: vetea, ranwez
"""

### Class representing a cluster of overlapping loci
class Cluster:
    __slots__ = ('name', 'loci_ref', 'loci_alt', 'end')
    def __init__(self, name, loci_ref=[], loci_alt=[], end=-1):
        self.name = name
        self.loci_ref = loci_ref
        self.loci_alt = loci_alt
        self.end = end
        
    def reverse_loci_coord(self):
        """
        Reverse the coordinates of all loci in the cluster.
        This is used for loci on the reverse strand.
        """
        for locus in self.loci_ref:
            locus.reverse(self.end)
        for locus in self.loci_alt:
            locus.reverse(self.end)

    def set_intervals(self) :
        """Pre compute the mRNA intervals for all loci in the cluster."""
        for locus in self.loci_ref:
            locus.set_mrna_intervals()
        for locus in self.loci_alt:
            locus.set_mrna_intervals()

        
        
        
        