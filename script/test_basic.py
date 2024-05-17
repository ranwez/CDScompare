#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:21:54 2024

@author: vetea
"""

import os, sys

# code adapted from https://csatlas.com/python-import-file-module/#import_a_file_in_a_different_directory
script_dir = os.path.dirname( __file__ )
sys.path.append( script_dir )
import main

## This script tests the 'main.py' program on multiple basic 'artificial' test files and checks if their return values match what is expected


# test function for the 'main.py' function 'get_gff_borders' (CDS coordinates acquisition function)
def test_get_gff_borders():
    
    # dictionary of inputs and expected ouputs for each test file for the 'main.py' function 'get_gff_borders' (CDS coordinates acquisition function)
    test_dict = {
        "basic" : ["./data/tests/basic_test.gff3", 
                        {'chr2A_00611930': [1,[100, 130, 150, 210, 240, 300]]}],
        
        "identical" : ["./data/tests/identical_test.gff3",
                            {'chr2A_00611930': [1,[100, 130, 150, 210, 240, 300]]}],
        
        "minus-CDS" : ["./data/tests/minus-CDS_test.gff3",
                       {'chr2A_00611930': [1,[100, 130, 240, 300]]}],
        
        "fusion" : ["./data/tests/fusion_test.gff3",
                    {'chr2A_00611930': [1,[100, 210, 240, 300]]}],
        
        "shift" : ["./data/tests/shift_test.gff3",
                   {'chr2A_00611930': [1,[100, 130, 151, 210, 240, 300]]}],
        
        "reverse" : ["./data/tests/reverse_test.gff3",
                     {'chr2A_00611930': [0,[300, 240, 210, 150, 130, 100]]}],
        
        "diff-start-before" : ["./data/tests/diff-start-before_test.gff3",
                               {'chr2A_00611930': [1,[40, 70, 90, 150, 180, 240]]}],
        
        "diff-start-after" : ["./data/tests/diff-start-after_test.gff3",
                              {'chr2A_00611930': [1,[160, 190, 210, 270, 300, 360]]}],
        
        "basic-2-loci" : ["./data/tests/basic-2-loci_test.gff3",
                              {'chr2A_00611930': [1,[100, 130, 150, 210, 240, 300]],
                               'chr2A_00620000': [1,[600, 700, 800, 900]]}],
        
        "identical-2-loci" : ["./data/tests/identical-2-loci_test.gff3",
                              {'chr2A_00611930': [1,[100, 130, 150, 210, 240, 300]],
                               'chr2A_00620000': [1,[600, 700, 800, 900]]}]
        }
    
    print("\n*************Testing the get_gff_borders function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        result = main.get_gff_borders(test_dict[test][0], False, False)
        print(result)
        print(test_dict[test][1])
        assert result == test_dict[test][1]
        
        
# test function for the 'main.py' function 'annotation_sort' (locus list creation and sorting function)
def test_annotation_sort():
    
    # dictionary of inputs and expected ouputs for each test file for the 'main.py' function 'annotation_sort' (locus list creation and sorting function)
    test_dict = {
        "overlapping-loci" : [{'chr_2A_1000' : [1,[50, 150]], 
                    'chr_2A_2000' : [1,[200, 350]],
                    'chr_2A_3000' : [1,[400, 550]],
                    'chr_2A_4000' : [1,[650, 700]],
                    'chr_2A_5000' : [1,[750, 800]]},
                              {'chr_2A_1000' : [1,[100, 250]],
                                          'chr_2A_2000' : [1,[300, 450]],
                                          'chr_2A_3000' : [1,[500, 600]],
                                          'chr_2A_4000' : [1,[650, 700]],
                                          'chr_2A_5000' : [1,[750, 780]],
                                          'chr_2A_6000' : [1,[790, 850]]},
                              [(50, 150, 'chr_2A_1000', True), 
                               (100, 250, 'chr_2A_1000', False), 
                               (200, 350, 'chr_2A_2000', True), 
                               (300, 450, 'chr_2A_2000', False), 
                               (400, 550, 'chr_2A_3000', True), 
                               (500, 600, 'chr_2A_3000', False), 
                               (650, 700, 'chr_2A_4000', False), 
                               (650, 700, 'chr_2A_4000', True), 
                               (750, 780, 'chr_2A_5000', False), 
                               (750, 800, 'chr_2A_5000', True), 
                               (790, 850, 'chr_2A_6000', False)]]
        }
    
    print("\n*************Testing the annotation_sort function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        result = main.annotation_sort(test_dict[test][0], test_dict[test][1], False, False)
        print(result)
        print(test_dict[test][2])
        assert result == test_dict[test][2]
        
        
# test function for the 'main.py' function 'locus_append_delete' (locus fusion in dictionary function)
def test_locus_append_delete():
    
    # dictionary of inputs and expected ouputs for each test file for the 'main.py' function 'locus_append_delete' (locus fusion in dictionary function)
    test_dict = {
        'simple' : [{'chr_2A_1000' : [1,[100,200, 300, 400]],
                    'chr_2A_2000' : [1,[500, 600, 700, 800]]},
                    'chr_2A_1000',
                    'chr_2A_2000',
                    {'chr_2A_1000' : [1,[100,200, 300, 400, 500, 600, 700, 800]]}]
        }
        
    print("\n*************Testing the locus_append_delete function*************")
    
    for test in test_dict:
        
        print(f"\n{test} test")
        
        main.locus_append_delete(test_dict[test][0], test_dict[test][1], test_dict[test][2], False, False)
        print(test_dict[test][0])
        print(test_dict[test][3])
        assert test_dict[test][0] == test_dict[test][3]
        
        
# test function for the 'main.py' function 'fuse_superloci' (overlapping loci fusion function)
def test_fuse_superloci():
    
    # dictionary of inputs and expected ouputs for each test file for the 'main.py' function 'fuse_superloci' (overlapping loci fusion function)
    test_dict = {
        'simple' : [{'chr_2A_1000' : [1,[50, 150]], 
                     'chr_2A_2000' : [1,[200, 350]],
                     'chr_2A_3000' : [1,[400, 550]],
                     'chr_2A_4000' : [1,[650, 700]],
                     'chr_2A_5000' : [1,[750, 800]]},
                    {'chr_2A_1000' : [1,[100, 250]],
                     'chr_2A_2000' : [1,[300, 450]],
                     'chr_2A_3000' : [1,[500, 600]],
                     'chr_2A_4000' : [1,[650, 700]],
                     'chr_2A_5000' : [1,[750, 780]],
                     'chr_2A_6000' : [1,[790, 850]]},
                    [(50, 150, 'chr_2A_1000', True), 
                     (100, 250, 'chr_2A_1000', False), 
                     (200, 350, 'chr_2A_2000', True), 
                     (300, 450, 'chr_2A_2000', False), 
                     (400, 550, 'chr_2A_3000', True), 
                     (500, 600, 'chr_2A_3000', False), 
                     (650, 700, 'chr_2A_4000', False), 
                     (650, 700, 'chr_2A_4000', True), 
                     (750, 780, 'chr_2A_5000', False), 
                     (750, 800, 'chr_2A_5000', True), 
                     (790, 850, 'chr_2A_6000', False)],
                    {'chr_2A_1000': [1, [50, 150, 200, 350, 400, 550]], 
                     'chr_2A_4000': [1, [650, 700]], 
                     'chr_2A_5000': [1, [750, 800]]},
                    {'chr_2A_1000': [1, [100, 250, 300, 450, 500, 600]], 
                     'chr_2A_4000': [1, [650, 700]], 
                     'chr_2A_5000': [1, [750, 780, 790, 850]]}]
        }

    print("\n*************Testing the fuse_superloci function*************")
    
    for test in test_dict:
        
        print(f"\n{test} test")
        
        main.fuse_superloci(test_dict[test][0], test_dict[test][1], test_dict[test][2], False, False)
        print(test_dict[test][0])
        print(test_dict[test][3])
        print("\n")
        print(test_dict[test][1])
        print(test_dict[test][4])
        assert test_dict[test][0] == test_dict[test][3]
        assert test_dict[test][1] == test_dict[test][4]
    

# test function for the 'main.py' function 'get_area_bounds' (CDS coordinates fusion function)
def test_get_area_bounds():
    
    # dictionary of inputs and expected ouputs for each test file for the 'main.py' function 'get_area_bounds' (CDS coordinates fusion function)
    test_dict = {
        "identical" : [[100, 130, 150, 210, 240, 300], 
                       [100, 130, 150, 210, 240, 300], 
                       [100, 130, 150, 210, 240, 300]],
        
        "minus-CDS" : [[100, 130, 150, 210, 240, 300],
                       [100, 130, 240, 300],
                       [100, 130, 150, 210, 240, 300]],
        
        "fusion" : [[100, 130, 150, 210, 240, 300],
                    [100, 210, 240, 300],
                    [100, 130, 150, 210, 240, 300]],
        
        "shift" : [[100, 130, 150, 210, 240, 300],
                   [100, 130, 151, 210, 240, 300],
                   [100, 130, 150, 151, 210, 240, 300]],
        
        "reverse" : [[300, 240, 210, 150, 130, 100],
                     [300, 240, 210, 150, 130, 100],
                     [300, 240, 210, 150, 130, 100]],
        
        "diff-start-before" : [[100, 130, 150, 210, 240, 300],
                               [40, 70, 90, 150, 180, 240],
                               [40, 70, 90, 100, 130, 150, 180, 210, 240, 300]],
        
        "diff-start-after" : [[100, 130, 150, 210, 240, 300],
                              [160, 190, 210, 270, 300, 360],
                              [100, 130, 150, 160, 190, 210, 240, 270, 300, 360]],
        
        "basic-2-loci (second locus)" : [[600, 700, 800, 900],
                                         [600, 700, 800, 900],
                                         [600, 700, 800, 900]]
        }
    
    print("\n*************Testing the get_area_bounds function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        result = main.get_area_bounds(test_dict[test][0], test_dict[test][1], False, False)
        print(result)
        print(test_dict[test][2])
        assert result == test_dict[test][2]
    

# test function for the 'main.py' function 'is_in_cds' (CDS inclusion in comparison area function)
def test_is_in_cds():
    
    # dictionary of inputs and expected ouputs for each test file for the 'main.py' function 'is_in_cds' (CDS inclusion in comparison area function)
    test_dict = {
        "identical" : [[100, 130, 150, 210, 240, 300],
                       [100, 130, 150, 210, 240, 300],
                       [True, False, True, False, True]],
        
        "minus-CDS" : [[100, 130, 240, 300],
                       [100, 130, 150, 210, 240, 300],
                       [True, False, False, False, True]],
        
        "fusion" : [[100, 210, 240, 300],
                    [100, 130, 150, 210, 240, 300],
                    [True, True, True, False, True]],
        
        "shift" : [[100, 130, 151, 210, 240, 300],
                   [100, 130, 150, 151, 210, 240, 300],
                   [True, False, False, True, False, True]],
        
        "reverse" : [[300, 240, 210, 150, 130, 100],
                     [300, 240, 210, 150, 130, 100],
                     [True, False, True, False, True]],
        
        "diff-start-before" : [[40, 70, 90, 150, 180, 240],
                               [40, 70, 90, 100, 130, 150, 180, 210, 240, 300],
                               [True, False, True, True, True, False, True, True, False]],
        
        "diff-start-after" : [[160, 190, 210, 270, 300, 360],
                              [100, 130, 150, 160, 190, 210, 240, 270, 300, 360],
                              [False, False, False, True, False, True, True, False, True]],
        
        "basic-2-loci (second locus)" : [[600, 700, 800, 900],
                                         [600, 700, 800, 900],
                                         [True, False, True]]
        }
    
    print("\n*************Testing the is_in_cds function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        result = main.is_in_cds(test_dict[test][0], test_dict[test][1], False, False)
        print(result)
        print(test_dict[test][2])
        assert result == test_dict[test][2]


# test function for the 'main.py' function 'compare_loci' (new locus comparison function)
def test_compare_loci():
    
    # dictionary of inputs and expected ouputs for each test file for the 'main.py' function 'compare_loci' (new locus comparison function)
    test_dict = {
        "identical" : [[100, 130, 150, 210, 240, 300], 
                       [100, 130, 150, 210, 240, 300],
                       [0, 150]],
        
        "minus-CDS" : [[100, 130, 150, 210, 240, 300],
                       [100, 130, 240, 300],
                       [60, 90]],
        
        "fusion" : [[100, 130, 150, 210, 240, 300],
                    [100, 210, 240, 300],
                    [140, 30]],
        
        "shift" : [[100, 130, 150, 210, 240, 300],
                   [100, 130, 151, 210, 240, 300],
                   [120, 30]],
        
        "reverse" : [[300, 240, 210, 150, 130, 100],
                     [300, 240, 210, 150, 130, 100],
                     [0, 150]],
        
        "diff-start-before" : [[100, 130, 150, 210, 240, 300],
                               [40, 70, 90, 150, 180, 240],
                               [210, 30]],
        
        "diff-start-after" : [[100, 130, 150, 210, 240, 300],
                              [160, 190, 210, 270, 300, 360],
                              [210, 30]],
        
        "basic-2-loci (second locus)" : [[600, 700, 800, 900],
                                         [600, 700, 800, 900],
                                         [0, 200]]
        }
    
    print("\n*************Testing the compare_loci function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        result = main.compare_loci(test_dict[test][0], test_dict[test][1], False, False)
        print(result)
        print(test_dict[test][2])
        assert result == test_dict[test][2]

    
# test function for the 'main.py' function 'create_vectors' (structure string creation function)
def test_create_vectors():
    
    # dictionary of inputs and expected ouputs for each test file for the 'main.py' function 'create_vectors' (structure string creation function)
    test_dict = {
        "basic" : [[100, 130, 150, 210, 240, 300], 
                   [100, "12312312312312312312312312312300000000000000000000123123123123123123123123123123123123123123123123123123123123000000000000000000000000000000123123123123123123123123123123123123123123123123123123123123"]],
        
        "identical" : [[100, 130, 150, 210, 240, 300], 
                       [100, "12312312312312312312312312312300000000000000000000123123123123123123123123123123123123123123123123123123123123000000000000000000000000000000123123123123123123123123123123123123123123123123123123123123"]],
        
        "minus-CDS" : [[100, 130, 240, 300], 
                       [100, "12312312312312312312312312312300000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000123123123123123123123123123123123123123123123123123123123123"]],
        
        "fusion" : [[100, 210, 240, 300], 
                    [100, "12312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312000000000000000000000000000000312312312312312312312312312312312312312312312312312312312312"]],
        
        "shift" : [[100, 130, 151, 210, 240, 300], 
                   [100, "12312312312312312312312312312300000000000000000000012312312312312312312312312312312312312312312312312312312312000000000000000000000000000000312312312312312312312312312312312312312312312312312312312312"]],
        
        "reverse" : [[300, 240, 210, 150, 130, 100], 
                     [100, "12312312312312312312312312312300000000000000000000123123123123123123123123123123123123123123123123123123123123000000000000000000000000000000123123123123123123123123123123123123123123123123123123123123"]],
        
        "diff-start-before" : [[40, 70, 90, 150, 180, 240], 
                               [40, "12312312312312312312312312312300000000000000000000123123123123123123123123123123123123123123123123123123123123000000000000000000000000000000123123123123123123123123123123123123123123123123123123123123"]],
        
        "diff-start-after" : [[160, 190, 210, 270, 300, 360],
                              [160, "12312312312312312312312312312300000000000000000000123123123123123123123123123123123123123123123123123123123123000000000000000000000000000000123123123123123123123123123123123123123123123123123123123123"]],
        
        "basic-2-loci (second locus)" : [[600, 700, 800, 900],
                                         [600, "123123123123123123123123123123123123123123123123123123123123123123123123123123123123123123123123123100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000002312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312"]]
        }
    
    print("\n*************Testing the create_vectors function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        result = main.create_vectors(test_dict[test][0], False, False)
        print(result)
        print(test_dict[test][1])
        assert result == test_dict[test][1]
     
        
# test function for the 'main.py' function 'create_vectors' (structure string creation function)
def test_old_compare_loci():
    
    # dictionary of inputs and expected ouputs for each test file for the 'main.py' function 'old_compare_loci' (old locus comparison function)
    test_dict = {
        "identical" : [[100, 130, 150, 210, 240, 300],
                       [100, 130, 150, 210, 240, 300],
                       [0, 150]],
        
        "minus-CDS" : [[100, 130, 150, 210, 240, 300],
                       [100, 130, 240, 300],
                       [60,90]],
        
        "fusion" : [[100, 130, 150, 210, 240, 300],
                    [100, 210, 240, 300],
                    [140,30]],
            
        "shift" : [[100, 130, 150, 210, 240, 300],
                   [100, 130, 151, 210, 240, 300],
                   [120,30]],
            
        "reverse" : [[100, 130, 150, 210, 240, 300],
                     [300, 240, 210, 150, 130, 100],
                     [0,150]],
            
        "diff-start-before" : [[100, 130, 150, 210, 240, 300],
                               [40, 70, 90, 150, 180, 240],
                               [210,30]],
            
        "diff-start-after" : [[100, 130, 150, 210, 240, 300],
                              [160, 190, 210, 270, 300, 360],
                              [210,30]],
        
        "basic-2-loci (second locus)" : [[600, 700, 800, 900],
                                         [600, 700, 800, 900],
                                         [0, 200]]
        }
    
    print("\n*************Testing the old_compare_loci function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        result = main.old_compare_loci(test_dict[test][0], test_dict[test][1], False, False)
        print(result)
        print(test_dict[test][2])
        assert result == test_dict[test][2]


# FUNCTIONS CALLS


# test function for the 'main.py' function 'get_gff_borders' (CDS coordinates acquisition function)
test_get_gff_borders()

# test function for the 'main.py' function 'get_area_bounds' (CDS coordinates fusion function)
test_get_area_bounds()

# test function for the 'main.py' function 'is_in_cds' (CDS inclusion in comparison area function)
test_is_in_cds()

# test function for the 'main.py' function 'compare_loci' (new locus comparison function)
test_compare_loci()
    
# test function for the 'main.py' function 'create_vectors' (structure string creation function)
test_create_vectors()
    
# test function for the 'main.py' function 'create_vectors' (structure string creation function)
test_old_compare_loci()

