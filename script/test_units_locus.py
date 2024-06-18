#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 10:59:29 2024

@author: vetea
"""

import os, sys

# code adapted from https://csatlas.com/python-import-file-module/#import_a_file_in_a_different_directory
script_dir = os.path.dirname( __file__ )
sys.path.append( script_dir )
import locus


# test function for the 'locus.py' class method 'reverse' (locus mRNAs inversion method)
def test_reverse():
    
    # dictionary of inputs and expected ouputs for the method
    test_dict = {
        
        "basic" : [[100, 129, 150, 209, 240, 299],
                   300,
                   {"test_mRNA": [1, 60, 91, 150, 171, 200]}],
        
        "diff_cluster_end" : [[100, 129, 150, 209, 240, 299],
                              400,
                              {"test_mRNA": [101, 160, 191, 250, 271, 300]}]
        
        }
    
    print("\n*************Testing the reverse function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        loc = locus.Locus(mRNAs = {"test_mRNA": test_dict[test][0]})
        loc.reverse(test_dict[test][1])
        result = loc.mRNAs
        print(f"result : {result}\n")
        print(test_dict[test][2])
        assert result == test_dict[test][2]