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
import annot_CSC
import locus
import cluster as cl
import intervals_utils as iu
import read_files as rf
import pre_comparison as pc
import comparison as comp

## This script tests the 'annot_CSC.py' program on multiple basic 'artificial' test files and checks if their return values match what is expected


# test function for the 'annot_CSC.py' function 'get_structure_id' (structure id acquisition function)
def test_get_structure_id():
    test_dict = {        
        "TRITD_HC" : [["chr2A", "PGSB_Jan2017", "gene", "22128", "24635", ".", "-", ".", "ID=TRITD2Av1G000030;primconf=HC\n"],
                      "TRITD2Av1G000030"],
        
        "EXP_chr2A" : [["chr2A", "exonerate:protein2genome:local", "gene", "608664", "611579", ".", "-", ".", "ID=chr2A_00611930;color=2;comment=modified:yes / Gene-Class:Non-canonical / start_changed_to3f;note=Origin:OSJnip_Chr04_04489784 / pred:prot2genome / prot-%25-ident:51.6 / prot-%25-cov:99.2118 / exo_corr:NA / Origin-Fam:NLR / Origin-Class:Canonical / noStart / Gene-Class:Non-canonical\n"],
                       "chr2A_00611930"],
        
        "chr2A" : [["chr2A", "exonerate:protein2genome:local", "gene", "608664", "611630", ".", "-", ".", "ID=chr2A_00611930;comment=Origin:OSJnip_Chr04_04489784 / pred:prot2genome / prot-%-ident:51.6 / prot-%-cov:99.2118 / exo_corr:NA / Origin-Fam:NLR / Origin-Class:Canonical / noStart / Gene-Class:Non-canonical;color=2\n"],
                   "chr2A_00611930"],
        
        "annot_best" : [["contig", "exonerate:protein2genome:local", "gene", "5131", "9328", ".", "+", ".", "ID=DWSvevo3_contig_0000005131;comment=Origin:DWSvevo1_chr2A_4125497 / pred:prot2genome / prot-%-ident:44.063 / prot-%-cov:68.2268 / score:51.7466 / scoreNC:46.9466 / exo_corr:corrected_intron / exo_corr:modif_stop / exo_corr:fix_overlap / Origin-Fam:Non-canonical / Origin-Class:LRR-RLK Gene-Class:Non-canonical / noStart / pbFrameshift / unexpectedSplicingSite / stopInFrame;color=2\n"],
                        "DWSvevo3_contig_0000005131"],
        
        "basic_test" : [["chr2A", "exonerate:protein2genome:local", "gene", "100", "300", ".", "-", ".", "ID=chr2A_00611930\n"],
                        "chr2A_00611930"]
        }
    
    print("\n*************Testing the get_structure_id function*************")
    
    for test in test_dict:
        
        print(f"\n{test} test")
        
        result = rf.get_structure_id(test_dict[test][0], False, False)
        print(f"result : {result}\n")
        print(test_dict[test][1])
        assert result == test_dict[test][1]


# test function for the 'annot_CSC.py' function 'get_gff_borders' (CDS coordinates acquisition function)
def test_get_gff_borders():
    
    # dictionary of inputs and expected ouputs for each test file for the 'annot_CSC.py' function 'get_gff_borders' (CDS coordinates acquisition function)
    test_dict = {
        "basic" : ["./data/tests/basic_test.gff3", 
                        {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}}],
        
        "identical" : ["./data/tests/identical_test.gff3",
                            {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}}],
        
        "minus-CDS" : ["./data/tests/minus-CDS_test.gff3",
                       {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 129, 240, 299]}}],
        
        "fusion" : ["./data/tests/fusion_test.gff3",
                    {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 209, 240, 299]}}],
        
        "shift" : ["./data/tests/shift_test.gff3",
                   {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 129, 151, 209, 240, 299]}}],
        
        "reverse" : ["./data/tests/reverse_test.gff3",
                     {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}}],
        
        "diff-start-before" : ["./data/tests/diff-start-before_test.gff3",
                               {"chr2A_00611930" : {'chr2A_00611930_mrna': [40, 69, 90, 149, 180, 239]}}],
        
        "diff-start-after" : ["./data/tests/diff-start-after_test.gff3",
                              {"chr2A_00611930" : {'chr2A_00611930_mrna': [160, 189, 210, 269, 300, 359]}}],
        
        "multiple mRNAs" : ["./data/tests/multiple_mRNAs_test.gff3",
                            {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 129, 150, 179, 240, 299],
                                                 'chr2A_00611930_mrna.2' : [100, 129, 240, 299]}}],
        
        "basic-2-loci" : ["./data/tests/basic-2-loci_test.gff3",
                              {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]},
                               "chr2A_00620000" : {'chr2A_00620000_mrna': [600, 699, 800, 899]}}],
        
        "identical-2-loci" : ["./data/tests/identical-2-loci_test.gff3",
                              {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]},
                               "chr2A_00620000" : {'chr2A_00620000_mrna': [600, 699, 800, 899]}}],
        
        "overlapping-loci" : ["./data/tests/overlapping-loci_test.gff3",
                              {"chr2A_1000" : {'chr2A_1000_mrna' : [50, 149]},
                               "chr2A_2000" : {'chr2A_2000_mrna' : [200, 349]},
                               "chr2A_3000" : {'chr2A_3000_mrna' : [400, 549]},
                              "chr2A_4000" : {'chr2A_4000_mrna' : [650, 699]},
                              "chr2A_5000" : {'chr2A_5000_mrna' : [750, 799]}}],
        
        "overlapping-loci-alt" : ["./data/tests/overlapping-loci-alt_test.gff3",
                                  {"chr2A_1000" : {'chr2A_1000_mrna' : [100, 249]},
                                   "chr2A_2000" : {'chr2A_2000_mrna' : [300, 449]},
                                   "chr2A_3000" : {'chr2A_3000_mrna' : [500, 599]},
                                   "chr2A_4000" : {'chr2A_4000_mrna' : [650, 699]},
                                   "chr2A_5000" : {'chr2A_5000_mrna' : [750, 779]},
                                   "chr2A_6000" : {'chr2A_6000_mrna' : [790, 849]}}],
        
        "length_computation_ref" : ["./data/tests/length_computation_ref_test.gff3", 
                        {"chr2A_00611930" : {'chr2A_00611930_mrna': [8, 13]}}],
        
        "length_computation_alt" : ["./data/tests/length_computation_alt_test.gff3", 
                        {"chr2A_00611930" : {'chr2A_00611930_mrna': [1, 4, 12, 13]}}]
        
        }
    
    print("\n*************Testing the get_gff_borders function*************")
    
    for test in test_dict:
        print(f"\n{test} file test")
        result = rf.get_gff_borders(test_dict[test][0], False, False)
        
        # verify presence of expected mRNAs
        for loc_id, loc in test_dict[test][1].items():
            print(f"Result keys : {result.keys()}")
            print(f"Expected key : {loc_id}")
            assert loc_id in result
            expected_mRNA = test_dict[test][1][loc_id]
            print(f"Result mRNAs : {result[loc_id].mRNAs}")
            print(f"Expected mRNA : {test_dict[test][1][loc_id]}")
            assert result[loc_id].contain_mrnas(**expected_mRNA)
        
        
# test function for the 'intervals_utils.py' class method 'transform_intervals_to_exclude_ub'
def test_transform_intervals_to_exclude_ub():
    
    # dictionary of inputs and expected ouputs for the method
    test_dict = {
        "basic" : [[100, 129, 150, 209, 240, 299], 
                   [100, 130, 150, 210, 240, 300]],
        
        "identical" : [[100, 129, 150, 209, 240, 299], 
                   [100, 130, 150, 210, 240, 300]],
        
        "minus-CDS" : [[100, 129, 240, 299], 
                   [100, 130, 240, 300]],
        
        "fusion" : [[100, 209, 240, 299], 
                   [100, 210, 240, 300]],
        
        "shift" : [[100, 129, 151, 209, 240, 299], 
                   [100, 130, 151, 210, 240, 300]],
        
        "diff-start-before" : [[40, 69, 90, 149, 180, 239], 
                               [40, 70, 90, 150, 180, 240]],
        
        "diff-start-after" : [[160, 189, 210, 269, 300, 359], 
                              [160, 190, 210, 270, 300, 360]],
        
        "basic-2-loci (second locus)" : [[600, 699, 800, 899], 
                                         [600, 700, 800, 900]],
        
        "length_computation_ref" : [[8, 13], 
                   [8, 14]],
        
        "length_computation_alt" : [[1, 4, 12, 13], 
                   [1, 5, 12, 14]]
        
        }
    
    print("\n*************Testing the transform_intervals_to_exclude_ub function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        # transform_intervals_to_exclude_ub is called by the class constructor
        interval = iu.OrderedIntervals(test_dict[test][0], True, False)
        result = interval.intervals
        print(f"result : {result}\n")
        print(test_dict[test][1])
        assert result == test_dict[test][1]
        
        
# test function for the 'intervals_utils.py' class method 'total_length'
def test_total_length():
    
    # dictionary of inputs and expected ouputs for the method
    test_dict = {
        "basic" : [[100, 129, 150, 209, 240, 299],
                   150],
        
        "identical" : [[100, 129, 150, 209, 240, 299],
                   150],
        
        "minus-CDS" : [[100, 129, 240, 299],
                   90],
        
        "fusion" : [[100, 209, 240, 299],
                   170],
        
        "shift" : [[100, 129, 151, 209, 240, 299],
                   149],
        
        "diff-start-before" : [[40, 69, 90, 149, 180, 239],
                               150],
        
        "diff-start-after" : [[160, 189, 210, 269, 300, 359],
                              150],
        
        "basic-2-loci (second locus)" : [[600, 699, 800, 899],
                                         200],
        
        "length_computation_ref" : [[8, 13],
                                    6],
        
        "length_computation_alt" : [[1, 4, 12, 13],
                                    6],
        
        }
    
    print("\n*************Testing the total_length function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        interval = iu.OrderedIntervals(test_dict[test][0], True, False)
        result = interval.total_length()
        print(f"result : {result}\n")
        print(test_dict[test][1])
        assert result == test_dict[test][1]
        
        
# test function for the 'intervals_utils.py' class method 'intersection'
def test_intersection():
    
    # dictionary of inputs and expected ouputs for the method
    test_dict = {
        
        "identical" : [[100, 129, 150, 209, 240, 299],
                       [100, 129, 150, 209, 240, 299],
                   [100, 129, 150, 209, 240, 299]],
        
        "minus-CDS" : [[100, 129, 150, 209, 240, 299],
                       [100, 129, 240, 299],
                   [100, 129, 240, 299]],
        
        "fusion" : [[100, 129, 150, 209, 240, 299],
                    [100, 209, 240, 299],
                   [100, 129, 150, 209, 240, 299]],
        
        "shift" : [[100, 129, 150, 209, 240, 299],
                   [100, 129, 151, 209, 240, 299],
                   [100, 129, 151, 209, 240, 299]],
        
        "diff-start-before" : [[100, 129, 150, 209, 240, 299],
                               [40, 69, 90, 149, 180, 239],
                               [100, 129, 180, 209]],
        
        "diff-start-after" : [[100, 129, 150, 209, 240, 299],
                              [160, 189, 210, 269, 300, 359],
                              [160, 189, 240, 269]],
        
        }
    
    print("\n*************Testing the intersection function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        interval1 = iu.OrderedIntervals(test_dict[test][0], True, False)
        interval2 = iu.OrderedIntervals(test_dict[test][1], True, False)
        result = interval1.intersection(interval2).get_intervals_with_included_ub()
        print(f"result : {result}\n")
        print(test_dict[test][2])
        assert result == test_dict[test][2]
        
        
# test function for the 'intervals_utils.py' class method 'union'
def test_union():
    
    # dictionary of inputs and expected ouputs for the method
    test_dict = {
        
        "identical" : [[100, 129, 150, 209, 240, 299],
                       [100, 129, 150, 209, 240, 299],
                   [100, 129, 150, 209, 240, 299]],
        
        "minus-CDS" : [[100, 129, 150, 209, 240, 299],
                       [100, 129, 240, 299],
                   [100, 129, 150, 209, 240, 299]],
        
        "fusion" : [[100, 129, 150, 209, 240, 299],
                    [100, 209, 240, 299],
                   [100, 209, 240, 299]],
        
        "shift" : [[100, 129, 150, 209, 240, 299],
                   [100, 129, 151, 209, 240, 299],
                   [100, 129, 150, 209, 240, 299]],
        
        "diff-start-before" : [[100, 129, 150, 209, 240, 299],
                               [40, 69, 90, 149, 180, 239],
                               [40, 69, 90, 299]],
        
        "diff-start-after" : [[100, 129, 150, 209, 240, 299],
                              [160, 189, 210, 269, 300, 359],
                              [100, 129, 150, 359]],
        
        }
    
    print("\n*************Testing the intersection function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        interval1 = iu.OrderedIntervals(test_dict[test][0], True, False)
        interval2 = iu.OrderedIntervals(test_dict[test][1], True, False)
        result = interval1.union(interval2).get_intervals_with_included_ub()
        print(f"result : {result}\n")
        print(test_dict[test][2])
        assert result == test_dict[test][2]
        
        
# test function for the 'locus.py' class method 'reverse'
def test_reverse():
    
    # dictionary of inputs and expected ouputs for the method
    test_dict = {
        
        "basic" : [[100, 129, 150, 209, 240, 299],
                   300,
                   {"test_mRNA": [1, 60, 91, 150, 171, 200]}],
        
        "diff_cluster" : [[100, 129, 150, 209, 240, 299],
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
        
        
# test function for the 'pre_comparison.py' class function 'get_reading_frame'
def test_get_reading_frame():
    
    # dictionary of inputs and expected ouputs for the method
    test_dict = {
        
        "identical" : [[100, 129, 150, 209, 240, 299],
                       [100, 129, 150, 209, 240, 299],
                       [1, 1, 1]],
        
        "minus-CDS_ref" : [[100, 129, 150, 209, 240, 299],
                           [100, 129, 240, 299],
                           [1, 1]],
        
        "minus-CDS_alt" : [[100, 129, 240, 299],
                           [100, 129, 240, 299],
                           [1, 1]],
        
        "fusion_ref" : [[100, 129, 150, 209, 240, 299],
                        [100, 129, 150, 209, 240, 299],
                        [1, 1, 1]],
        
        "fusion_alt" : [[100, 209, 240, 299],
                        [100, 129, 150, 209, 240, 299],
                        [1, 3, 3]],
        
        "shift_ref" : [[100, 129, 150, 209, 240, 299],
                       [100, 129, 151, 209, 240, 299],
                       [1, 2, 1]],
        
        "shift_alt" : [[100, 129, 151, 209, 240, 299],
                       [100, 129, 151, 209, 240, 299],
                       [1, 1, 3]],
        
        "diff-start-before_ref" : [[100, 129, 150, 209, 240, 299],
                                   [100, 129, 180, 209],
                                   [1, 1]],
        
        "diff-start-before_alt" : [[40, 69, 90, 149, 180, 239],
                                   [100, 129, 180, 209],
                                   [2, 1]],
        
        "diff-start-after_ref" : [[100, 129, 150, 209, 240, 299],
                                  [160, 189, 240, 269],
                                  [2, 1]],
        
        "diff-start-after_alt" : [[160, 189, 210, 269, 300, 359],
                                  [160, 189, 240, 269],
                                  [1, 1]],
        
        }
    
    print("\n*************Testing the get_reading_frame function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        result = pc.get_reading_frame(test_dict[test][0], test_dict[test][1])
        print(f"result : {result}\n")
        print(test_dict[test][2])
        assert result == test_dict[test][2]
        
        
# test function for the 'annot_CSC.py' function 'annotation_sort' (locus list creation and sorting function)
def test_annotation_sort():
    
    # dictionary of inputs and expected ouputs for each test file for the 'annot_CSC.py' function 'annotation_sort' (locus list creation and sorting function)
    test_dict = {"overlapping-loci" : [
                                {'chr2A_1000': locus.Locus(name='chr2A_1000', 
                                                     mRNAs={'chr2A_1000_mrna': [50, 149]}, 
                                                     start=50, 
                                                     end=149, 
                                                     direction='direct'), 
                                'chr2A_2000': locus.Locus(name='chr2A_2000', 
                                                     mRNAs={'chr2A_2000_mrna': [200, 349]}, 
                                                     start=200, 
                                                     end=349, 
                                                     direction='direct'), 
                                'chr2A_3000': locus.Locus(name='chr2A_3000', 
                                                     mRNAs={'chr2A_3000_mrna': [400, 549]}, 
                                                     start=400, 
                                                     end=549, 
                                                     direction='direct'),
                                'chr2A_4000': locus.Locus(name='chr2A_4000', 
                                                     mRNAs={'chr2A_4000_mrna': [650, 699]}, 
                                                     start=650, 
                                                     end=699, 
                                                     direction='direct'),
                                'chr2A_5000': locus.Locus(name='chr2A_5000', 
                                                     mRNAs={'chr2A_5000_mrna': [750, 799]}, 
                                                     start=750, 
                                                     end=799, 
                                                     direction='direct')},
                                                    
                                {'chr2A_1000': locus.Locus(name='chr2A_1000', 
                                                     mRNAs={'chr2A_1000_mrna': [100, 249]}, 
                                                     start=100, 
                                                     end=249, 
                                                     direction='direct'), 
                                'chr2A_2000': locus.Locus(name='chr2A_2000', 
                                                     mRNAs={'chr2A_2000_mrna': [300, 449]}, 
                                                     start=300, 
                                                     end=449, 
                                                     direction='direct'), 
                                'chr2A_3000': locus.Locus(name='chr2A_3000', 
                                                     mRNAs={'chr2A_3000_mrna': [500, 599]}, 
                                                     start=500, 
                                                     end=599, 
                                                     direction='direct'),
                                'chr2A_4000': locus.Locus(name='chr2A_4000', 
                                                     mRNAs={'chr2A_4000_mrna': [650, 699]}, 
                                                     start=650, 
                                                     end=699, 
                                                     direction='direct'),
                                'chr2A_5000': locus.Locus(name='chr2A_5000', 
                                                     mRNAs={'chr2A_5000_mrna': [750, 779]}, 
                                                     start=750, 
                                                     end=779, 
                                                     direction='direct'),
                                'chr2A_6000': locus.Locus(name='chr2A_6000', 
                                                     mRNAs={'chr2A_6000_mrna': [790, 849]}, 
                                                     start=790, 
                                                     end=849, 
                                                     direction='direct')},
                                                      
                                [(50, 149, 'chr2A_1000', True,'direct'), 
                               (100, 249, 'chr2A_1000', False,'direct'), 
                               (200, 349, 'chr2A_2000', True,'direct'), 
                               (300, 449, 'chr2A_2000', False,'direct'), 
                               (400, 549, 'chr2A_3000', True,'direct'), 
                               (500, 599, 'chr2A_3000', False,'direct'), 
                               (650, 699, 'chr2A_4000', False,'direct'), 
                               (650, 699, 'chr2A_4000', True,'direct'), 
                               (750, 779, 'chr2A_5000', False,'direct'), 
                               (750, 799, 'chr2A_5000', True,'direct'), 
                               (790, 849, 'chr2A_6000', False,'direct')]]
    }
    
    print("\n*************Testing the annotation_sort function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        result = pc.annotation_sort(test_dict[test][0], test_dict[test][1], False, False)
        print(f"result : {result}\n")
        print(f"expected result : {test_dict[test][2]}\n")
        assert result == test_dict[test][2]
        
        
# test function for the 'annot_CSC.py' function 'construct_clusters' (overlapping loci grouping function)
def test_construct_clusters():
    
    # dictionary of inputs and expected ouputs for each test file for the 'annot_CSC.py' function 'construct_clusters' (overlapping loci grouping function)
    test_dict = {'simple' : [{
        'chr2A_1000': locus.Locus(name='chr2A_1000', 
                                 mRNAs={'chr2A_1000_mrna': [50, 149]}, 
                                 start=50, 
                                 end=149, 
                                 direction='direct'), 
        'chr2A_2000': locus.Locus(name='chr2A_2000', 
                                 mRNAs={'chr2A_2000_mrna': [200, 349]}, 
                                 start=200, 
                                 end=349, 
                                 direction='direct'), 
        'chr2A_3000': locus.Locus(name='chr2A_3000', 
                                 mRNAs={'chr2A_3000_mrna': [400, 549]}, 
                                 start=400, 
                                 end=549, 
                                 direction='direct'),
        'chr2A_4000': locus.Locus(name='chr2A_4000', 
                                 mRNAs={'chr2A_4000_mrna': [650, 699]}, 
                                 start=650, 
                                 end=699, 
                                 direction='direct'),
        'chr2A_5000': locus.Locus(name='chr2A_5000', 
                                 mRNAs={'chr2A_5000_mrna': [750, 799]}, 
                                 start=750, 
                                 end=799, 
                                 direction='direct')},
                    
        {'chr2A_1000': locus.Locus(name='chr2A_1000', 
                             mRNAs={'chr2A_1000_mrna': [100, 249]}, 
                             start=100, 
                             end=249, 
                             direction='direct'), 
        'chr2A_2000': locus.Locus(name='chr2A_2000', 
                             mRNAs={'chr2A_2000_mrna': [300, 449]}, 
                             start=300, 
                             end=449, 
                             direction='direct'), 
        'chr2A_3000': locus.Locus(name='chr2A_3000', 
                             mRNAs={'chr2A_3000_mrna': [500, 599]}, 
                             start=500, 
                             end=599, 
                             direction='direct'),
        'chr2A_4000': locus.Locus(name='chr2A_4000', 
                             mRNAs={'chr2A_4000_mrna': [650, 699]}, 
                             start=650, 
                             end=699, 
                             direction='direct'),
        'chr2A_5000': locus.Locus(name='chr2A_5000', 
                             mRNAs={'chr2A_5000_mrna': [750, 779]}, 
                             start=750, 
                             end=779, 
                             direction='direct'),
        'chr2A_6000': locus.Locus(name='chr2A_6000', 
                             mRNAs={'chr2A_6000_mrna': [790, 849]}, 
                             start=790, 
                             end=849, 
                             direction='direct')},
        
        [(50, 149, 'chr2A_1000', True,'direct'), 
       (100, 249, 'chr2A_1000', False,'direct'), 
       (200, 349, 'chr2A_2000', True,'direct'), 
       (300, 449, 'chr2A_2000', False,'direct'), 
       (400, 549, 'chr2A_3000', True,'direct'), 
       (500, 599, 'chr2A_3000', False,'direct'), 
       (650, 699, 'chr2A_4000', False,'direct'), 
       (650, 699, 'chr2A_4000', True,'direct'), 
       (750, 779, 'chr2A_5000', False,'direct'), 
       (750, 799, 'chr2A_5000', True,'direct'), 
       (790, 849, 'chr2A_6000', False,'direct')],
        
        {'cluster 0': {'ref': [{'chr2A_1000_mrna': [50, 149]}, {'chr2A_2000_mrna': [200, 349]}, {'chr2A_3000_mrna': [400, 549]}], 'alt': [{'chr2A_1000_mrna': [100, 249]}, {'chr2A_2000_mrna': [300, 449]}, {'chr2A_3000_mrna': [500, 599]}]},
         'cluster 6': {'ref': [{'chr2A_4000_mrna': [650, 699]}], 'alt': [{'chr2A_4000_mrna': [650, 699]}]},
         'cluster 8': {'ref': [{'chr2A_5000_mrna': [750, 799]}], 'alt': [{'chr2A_5000_mrna': [750, 779]}, {'chr2A_6000_mrna': [790, 849]}]}}],
        
                "reverse-basic" : [{"chr2A_00611930" : locus.Locus(name='chr2A_00611930', 
                                     mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                     start=100, 
                                     end=299, 
                                     direction='direct')},
            
                                   {"chr2A_00611930" : locus.Locus(name='chr2A_00611930', 
                                     mRNAs={'chr2A_00611930_mrna': [299, 240, 209, 150, 129, 100]}, 
                                     start=299, 
                                     end=100, 
                                     direction='reverse')},
                                   
                                   [(100, 299, 'chr2A_00611930', True,'direct'), 
                                  (299, 100, 'chr2A_00611930', False,'reverse')], 
            
                                   {'cluster 0':  {'ref': [{'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}], 'alt': []},
                                    'cluster 1':  {'ref': [], 'alt': [{'chr2A_00611930_mrna': [299, 240, 209, 150, 129, 100]}]}}]
        
    }

    print("\n*************Testing the construct_clusters function*************")
    
    for test in test_dict:
        
        print(f"\n{test} test")
        
        result = pc.construct_clusters(test_dict[test][0], test_dict[test][1], test_dict[test][2], False, False)
        for cluster_id, cluster in test_dict[test][3].items():
            result_details = result[cluster_id].get_details()
            print(f"expected mRNAs: {test_dict[test][3][cluster_id]}")
            print(f"result mRNAs: {result_details}")
            assert result_details == cluster


# test function for the 'annot_CSC.py' function 'compare_loci' (new locus comparison function)
def test_compare_loci():
    
    # dictionary of inputs and expected ouputs for each test file for the 'annot_CSC.py' function 'compare_loci' (new locus comparison function)
    test_dict = {
        "identical" : [locus.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                 start=100, 
                                 end=299, 
                                 direction='direct'), 
                       locus.Locus(name='chr2A_00611930', 
                                    mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                    start=100, 
                                    end=299, 
                                    direction='direct'),
                       ([150, 0, 0], 100.0, ([], [[]]), 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
        "minus-CDS" : [locus.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                 start=100, 
                                 end=299, 
                                 direction='direct'),
                       locus.Locus(name='chr2A_00611930', 
                                    mRNAs={'chr2A_00611930_mrna': [100, 129, 240, 299]}, 
                                    start=100, 
                                    end=299, 
                                    direction='direct'),
                       ([90, 60, 0], 60.0, ([150, 210], [[]]), 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
        "fusion" : [locus.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                 start=100, 
                                 end=299, 
                                 direction='direct'),
                    locus.Locus(name='chr2A_00611930', 
                                mRNAs={'chr2A_00611930_mrna': [100, 209, 240, 299]}, 
                                start=100, 
                                end=299, 
                                direction='direct'),
                    ([30, 20, 120], 17.6, ([130, 150], [[150, 209], [240, 299]]), 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
        "shift" : [locus.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                 start=100, 
                                 end=299, 
                                 direction='direct'),
                   locus.Locus(name='chr2A_00611930', 
                                mRNAs={'chr2A_00611930_mrna': [100, 129, 151, 209, 240, 299]}, 
                                start=100, 
                                end=299, 
                                direction='direct'),
                   ([30, 1, 119], 20.0, ([150, 151], [[151, 209], [240, 299]]), 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
        "reverse-basic" : [locus.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                 start=100, 
                                 end=299, 
                                 direction='direct'),
                     locus.Locus(name='chr2A_00611930', 
                                mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                start=100, 
                                end=299, 
                                direction='reverse'),
                     ('_', 0.0, '_', '_', '_')],
        
        "diff-start-before" : [locus.Locus(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                            start=100, 
                                            end=299, 
                                            direction='direct'),
                               locus.Locus(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [40, 69, 90, 149, 180, 239]}, 
                                            start=40, 
                                            end=239, 
                                            direction='direct'),
                               ([30, 180, 30], 12.5, ([40, 70, 90, 100, 130, 180, 210, 300], [[100, 129]]), 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
        "diff-start-after" : [locus.Locus(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                            start=100, 
                                            end=299, 
                                            direction='direct'),
                              locus.Locus(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [160, 189, 210, 269, 300, 359]}, 
                                            start=160, 
                                            end=359, 
                                            direction='direct'),
                              ([30, 180, 30], 12.5, ([100, 130, 150, 160, 190, 240, 270, 360], [[160, 189]]), 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
        "multiple_mRNAs" : [locus.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                 start=100, 
                                 end=299, 
                                 direction='direct'), 
                       locus.Locus(name='chr2A_00611930', 
                                    mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 179, 240, 299],
                                           'chr2A_00611930_mrna.2': [100, 129, 240, 299]}, 
                                    start=100, 
                                    end=299, 
                                    direction='direct'),
                       ([120, 30, 0], 80.0, ([180, 210], [[]]), 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
        "basic-2-loci (second locus)" : [locus.Locus(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [600, 699, 800, 899]}, 
                                            start=600, 
                                            end=899, 
                                            direction='direct'),
                                         locus.Locus(name='chr2A_00611930', 
                                                    mRNAs={'chr2A_00611930_mrna': [600, 699, 800, 899]}, 
                                                    start=600, 
                                                    end=899, 
                                                    direction='direct'),
                                         ([200, 0, 0], 100.0, ([], [[]]), 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
        "length_computation" : [locus.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [8, 13]}, 
                                 start=8, 
                                 end=13, 
                                 direction='direct'), 
                       locus.Locus(name='chr2A_00611930', 
                                    mRNAs={'chr2A_00611930_mrna': [1, 4, 12, 13]}, 
                                    start=1, 
                                    end=13, 
                                    direction='direct'),
                       ([2, 8, 0], 20.0, ([1, 5, 8, 12], [[]]), 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')]
        
        }
    
    print("\n*************Testing the compare_loci function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        result = comp.compare_loci(test_dict[test][0], test_dict[test][1], False, False)
        print(f"Expected result = {test_dict[test][2]}")
        print(f"Result = {result}")
        assert result == test_dict[test][2]

    
# test function for the 'annot_CSC.py' function 'create_vectors' (structure string creation function)
def test_create_vectors():
    
    # dictionary of inputs and expected ouputs for each test file for the 'annot_CSC.py' function 'create_vectors' (structure string creation function)
    test_dict = {
        "basic" : [[100, 129, 150, 209, 240, 299], 
                   [100, "12312312312312312312312312312300000000000000000000123123123123123123123123123123123123123123123123123123123123000000000000000000000000000000123123123123123123123123123123123123123123123123123123123123"]],
        
        "identical" : [[100, 129, 150, 209, 240, 299], 
                       [100, "12312312312312312312312312312300000000000000000000123123123123123123123123123123123123123123123123123123123123000000000000000000000000000000123123123123123123123123123123123123123123123123123123123123"]],
        
        "minus-CDS" : [[100, 129, 240, 299], 
                       [100, "12312312312312312312312312312300000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000123123123123123123123123123123123123123123123123123123123123"]],
        
        "fusion" : [[100, 209, 240, 299], 
                    [100, "12312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312000000000000000000000000000000312312312312312312312312312312312312312312312312312312312312"]],
        
        "shift" : [[100, 129, 151, 209, 240, 299], 
                   [100, "12312312312312312312312312312300000000000000000000012312312312312312312312312312312312312312312312312312312312000000000000000000000000000000312312312312312312312312312312312312312312312312312312312312"]],
        
        "diff-start-before" : [[40, 69, 90, 149, 180, 239], 
                               [40, "12312312312312312312312312312300000000000000000000123123123123123123123123123123123123123123123123123123123123000000000000000000000000000000123123123123123123123123123123123123123123123123123123123123"]],
        
        "diff-start-after" : [[160, 189, 210, 269, 300, 359],
                              [160, "12312312312312312312312312312300000000000000000000123123123123123123123123123123123123123123123123123123123123000000000000000000000000000000123123123123123123123123123123123123123123123123123123123123"]],
        
        "basic-2-loci (second locus)" : [[600, 699, 800, 899],
                                         [600, "123123123123123123123123123123123123123123123123123123123123123123123123123123123123123123123123123100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000002312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312"]],
        
        "length_computation_ref" : [[8, 13], 
                   [8, "123123"]],
        
        "length_computation_alt" : [[1, 4, 12, 13], 
                   [1, "1231000000023"]]
        
        }
    
    print("\n*************Testing the create_vectors function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        result = pc.create_vectors(test_dict[test][0], False, False)
        print(f"result : {result}\n")
        print(test_dict[test][1])
        assert result == test_dict[test][1]
     
        
# test function for the 'annot_CSC.py' function 'create_vectors' (structure string creation function)
def test_old_compare_loci():
    
    # dictionary of inputs and expected ouputs for each test file for the 'annot_CSC.py' function 'old_compare_loci' (old locus comparison function)
    test_dict = {
        "identical" : [locus.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                 start=100, 
                                 end=299, 
                                 direction='direct'), 
                       locus.Locus(name='chr2A_00611930', 
                                    mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                    start=100, 
                                    end=299, 
                                    direction='direct'),
                       ([150, 0, 0], 100.0, 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
        "minus-CDS" : [locus.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                 start=100, 
                                 end=299, 
                                 direction='direct'),
                       locus.Locus(name='chr2A_00611930', 
                                    mRNAs={'chr2A_00611930_mrna': [100, 129, 240, 299]}, 
                                    start=100, 
                                    end=299, 
                                    direction='direct'),
                       ([90, 60, 0], 60.0, 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
        "fusion" : [locus.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                 start=100, 
                                 end=299, 
                                 direction='direct'),
                    locus.Locus(name='chr2A_00611930', 
                                mRNAs={'chr2A_00611930_mrna': [100, 209, 240, 299]}, 
                                start=100, 
                                end=299, 
                                direction='direct'),
                    ([30, 20, 120], 17.6, 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
        "shift" : [locus.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                 start=100, 
                                 end=299, 
                                 direction='direct'),
                   locus.Locus(name='chr2A_00611930', 
                                mRNAs={'chr2A_00611930_mrna': [100, 129, 151, 209, 240, 299]}, 
                                start=100, 
                                end=299, 
                                direction='direct'),
                   ([30, 1, 119], 20.0, 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],

        "reverse-basic" : [locus.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                 start=100, 
                                 end=299, 
                                 direction='direct'),
                     locus.Locus(name='chr2A_00611930', 
                                mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                start=100, 
                                end=299, 
                                direction='reverse'),
                     ('_', 0.0, '_', '_')],
        
        "reverse-reverse" : [locus.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                 start=100, 
                                 end=299, 
                                 direction='reverse'),
                     locus.Locus(name='chr2A_00611930', 
                                mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                start=100, 
                                end=299, 
                                direction='reverse'),
                     ([150, 0, 0], 100.0, 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
        "diff-start-before" : [locus.Locus(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                            start=100, 
                                            end=299, 
                                            direction='direct'),
                               locus.Locus(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [40, 69, 90, 149, 180, 239]}, 
                                            start=40, 
                                            end=239, 
                                            direction='direct'),
                               ([30, 180, 30], 12.5, 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
        "diff-start-after" : [locus.Locus(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                            start=100, 
                                            end=299, 
                                            direction='direct'),
                              locus.Locus(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [160, 189, 210, 269, 300, 359]}, 
                                            start=160, 
                                            end=359, 
                                            direction='direct'),
                              ([30, 180, 30], 12.5, 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
        "multiple_mRNAs" : [locus.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                 start=100, 
                                 end=299, 
                                 direction='direct'), 
                       locus.Locus(name='chr2A_00611930', 
                                    mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 179, 240, 299],
                                           'chr2A_00611930_mrna.2': [100, 129, 240, 299]}, 
                                    start=100, 
                                    end=299, 
                                    direction='direct'),
                       ([120, 30, 0], 80.0, 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
        "basic-2-loci (second locus)" : [locus.Locus(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [600, 699, 800, 899]}, 
                                            start=600, 
                                            end=899, 
                                            direction='direct'),
                                         locus.Locus(name='chr2A_00611930', 
                                                    mRNAs={'chr2A_00611930_mrna': [600, 699, 800, 899]}, 
                                                    start=600, 
                                                    end=899, 
                                                    direction='direct'),
                                         ([200, 0, 0], 100.0, 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
        'overlapping-loci (first cluster/locus)' : [locus.Locus(name='chr2A_1000', 
                                                             mRNAs={'chr2A_1000_mrna': [50, 149]}, 
                                                             start=50, 
                                                             end=149, 
                                                             direction='direct'),
                                                    locus.Locus(name='chr2A_1000', 
                                                              mRNAs={'chr2A_1000_mrna': [100, 249]}, 
                                                              start=100, 
                                                              end=249, 
                                                              direction='direct'),
                                                    ([0, 150, 50], 0.0, 'chr2A_1000_mrna', 'chr2A_1000_mrna')],
        
        'overlapping-loci (first cluster / 2nd locus)' : [locus.Locus(name='chr2A_2000', 
                                                             mRNAs={'chr2A_2000_mrna': [200, 349]}, 
                                                             start=200, 
                                                             end=349, 
                                                             direction='direct'),
                                                          locus.Locus(name='chr2A_2000', 
                                                              mRNAs={'chr2A_2000_mrna': [300, 449]}, 
                                                              start=300, 
                                                              end=449, 
                                                              direction='direct'),
                                                          ([0, 200, 50], 0.0, 'chr2A_2000_mrna', 'chr2A_2000_mrna')],
                         
        "length_computation" : [locus.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [8, 13]}, 
                                 start=8, 
                                 end=13, 
                                 direction='direct'), 
                       locus.Locus(name='chr2A_00611930', 
                                    mRNAs={'chr2A_00611930_mrna': [1, 4, 12, 13]}, 
                                    start=1, 
                                    end=13, 
                                    direction='direct'),
                       ([2, 8, 0], 20.0, 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')]
        }
    
    print("\n*************Testing the old_compare_loci function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        result = comp.old_compare_loci(test_dict[test][0], test_dict[test][1], False, False)
        print(f"result : {result}\n")
        print(f"expected result : {test_dict[test][2]}\n")
        assert result == test_dict[test][2]

        
# test function for the 'annot_CSC.py' function 'annotation_match' (main annotation comparison function) with the 'create_strings' parameter as 'False' (uses new program version)
def test_new_annotation_match():
    
    # dictionary of inputs and expected ouputs for each test file for the 'annot_CSC.py' function 'annotation_match' (main annotation comparison function)
    test_dict = {'basic' : [cl.Cluster(name="cluster 0",
                                       loci={'ref': [locus.Locus(name='chr2A_00611930', 
                                                                    mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                                                    start=100, 
                                                                    end=299, 
                                                                    direction='direct')],
                                             'alt': [locus.Locus(name='chr2A_00611930', 
                                                                          mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                                                          start=100, 
                                                                          end=299, 
                                                                          direction='direct')]}),
                             [{"reference" : 'chr2A_00611930',
                                "reference start" : 100,
                                "reference end" : 299,
                                "alternative" : 'chr2A_00611930',
                                "alternative start" : 100,
                                "alternative end" : 299,
                                "reference mRNA" : 'chr2A_00611930_mrna',
                                "alternative mRNA" : 'chr2A_00611930_mrna',
                                "mismatch/match" : [150, 0, 0],
                                "identity" : 100.0,
                                "mismatch zones" : ([], [[]]),
                                "cluster name" : "cluster 0",
                                "reference mRNA number" : 1,
                                "alternative mRNA number" : 1}]],
                 
                 'shift' : [cl.Cluster(name="cluster 0",
                                                    loci={'ref': [locus.Locus(name='chr2A_00611930', 
                                                                                 mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                                                                 start=100, 
                                                                                 end=299, 
                                                                                 direction='direct')],
                                                          'alt': [locus.Locus(name='chr2A_00611930', 
                                                                                       mRNAs={'chr2A_00611930_mrna': [100, 129, 151, 209, 240, 299]}, 
                                                                                       start=100, 
                                                                                       end=299, 
                                                                                       direction='direct')]}),
                            [{"reference" : 'chr2A_00611930',
                            "reference start" : 100,
                            "reference end" : 299,
                            "alternative" : 'chr2A_00611930',
                            "alternative start" : 100,
                            "alternative end" : 299,
                            "reference mRNA" : 'chr2A_00611930_mrna',
                            "alternative mRNA" : 'chr2A_00611930_mrna',
                            "mismatch/match" : [30, 1, 119],
                            "identity" : 20.0,
                            "mismatch zones" : ([150, 151], [[151, 209], [240, 299]]),
                            "cluster name" : "cluster 0",
                            "reference mRNA number" : 1,
                            "alternative mRNA number" : 1}]],
                 
                 'multiple_mRNAs' : [cl.Cluster(name="cluster 0",
                                                loci={'ref': [locus.Locus(name='chr2A_00611930', 
                                                                         mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                                                         start=100, 
                                                                         end=299, 
                                                                         direction='direct')],
                                                      'alt': [locus.Locus(name='chr2A_00611930', 
                                                                         mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 179, 240, 299],
                                                                                'chr2A_00611930_mrna.2': [100, 129, 240, 299]}, 
                                                                         start=100, 
                                                                         end=299, 
                                                                         direction='direct')]}),
                                          [{"reference" : 'chr2A_00611930',
                                             "reference start" : 100,
                                             "reference end" : 299,
                                             "alternative" : 'chr2A_00611930',
                                             "alternative start" : 100,
                                             "alternative end" : 299,
                                             "reference mRNA" : 'chr2A_00611930_mrna',
                                             "alternative mRNA" : 'chr2A_00611930_mrna',
                                             "mismatch/match" : [120, 30, 0],
                                             "identity" : 80.0,
                                             "mismatch zones" : ([180, 210], [[]]),
                                             "cluster name" : "cluster 0",
                                             "reference mRNA number" : 1,
                                             "alternative mRNA number" : 2}]],
                 
                 'overlapping-loci (first cluster)' : [cl.Cluster(name="cluster 0",
                                                                  loci={'ref': [locus.Locus(name='chr2A_1000', 
                                              mRNAs={'chr2A_1000_mrna': [50, 149]}, 
                                              start=50, 
                                              end=149, 
                                              direction='direct'), 
                                        locus.Locus(name='chr2A_2000', 
                                              mRNAs={'chr2A_2000_mrna': [200, 349]}, 
                                              start=200, 
                                              end=349, 
                                              direction='direct'), 
                                        locus.Locus(name='chr2A_3000', 
                                              mRNAs={'chr2A_3000_mrna': [400, 549]}, 
                                              start=400, 
                                              end=549, 
                                              direction='direct')],
                                                                        'alt': [locus.Locus(name='chr2A_1000', 
                                          mRNAs={'chr2A_1000_mrna': [100, 249]}, 
                                          start=100, 
                                          end=249, 
                                          direction='direct'), 
                                        locus.Locus(name='chr2A_2000', 
                                          mRNAs={'chr2A_2000_mrna': [300, 449]}, 
                                          start=300, 
                                          end=449, 
                                          direction='direct'), 
                                        locus.Locus(name='chr2A_3000', 
                                          mRNAs={'chr2A_3000_mrna': [500, 599]}, 
                                          start=500, 
                                          end=599, 
                                          direction='direct')]}),
                                    [{"reference" : 'chr2A_1000',
                                   "reference start" : 50,
                                   "reference end" : 149,
                                   "alternative" : 'chr2A_1000',
                                   "alternative start" : 100,
                                   "alternative end" : 249,
                                   "reference mRNA" : 'chr2A_1000_mrna',
                                   "alternative mRNA" : 'chr2A_1000_mrna',
                                   "mismatch/match" : [0, 150, 50],
                                   "identity" : 0.0,
                                   "mismatch zones" : ([50, 100, 150, 250], [[100, 149]]),
                                   "cluster name" : "cluster 0",
                                   "reference mRNA number" : 1,
                                   "alternative mRNA number" : 1},
                                     {"reference" : 'chr2A_2000',
                                    "reference start" : 200,
                                    "reference end" : 349,
                                    "alternative" : 'chr2A_2000',
                                    "alternative start" : 300,
                                    "alternative end" : 449,
                                    "reference mRNA" : 'chr2A_2000_mrna',
                                    "alternative mRNA" : 'chr2A_2000_mrna',
                                    "mismatch/match" : [0, 200, 50],
                                    "identity" : 0.0,
                                    "mismatch zones" : ([200, 300, 350, 450], [[300, 349]]),
                                    "cluster name" : "cluster 0",
                                    "reference mRNA number" : 1,
                                    "alternative mRNA number" : 1},
                                     {"reference" : 'chr2A_3000',
                                     "reference start" : 400,
                                     "reference end" : 549,
                                     "alternative" : 'chr2A_3000',
                                     "alternative start" : 500,
                                     "alternative end" : 599,
                                     "reference mRNA" : 'chr2A_3000_mrna',
                                     "alternative mRNA" : 'chr2A_3000_mrna',
                                     "mismatch/match" : [0, 150, 50],
                                     "identity" : 0.0,
                                     "mismatch zones" : ([400, 500, 550, 600], [[500, 549]]),
                                     "cluster name" : "cluster 0",
                                     "reference mRNA number" : 1,
                                     "alternative mRNA number" : 1}]],
                 
                 'length_computation' : [cl.Cluster(name="cluster 0",
                                                    loci={'ref': [locus.Locus(name='chr2A_00611930', 
                                          mRNAs={'chr2A_00611930_mrna': [8, 13]}, 
                                          start=8, 
                                          end=13, 
                                          direction='direct')], 
                                                          'alt': [locus.Locus(name='chr2A_00611930', 
                                             mRNAs={'chr2A_00611930_mrna': [1, 4, 12, 13]}, 
                                             start=1, 
                                             end=13, 
                                             direction='direct')]}),
                                          [{"reference" : 'chr2A_00611930',
                                             "reference start" : 8,
                                             "reference end" : 13,
                                             "alternative" : 'chr2A_00611930',
                                             "alternative start" : 1,
                                             "alternative end" : 13,
                                             "reference mRNA" : 'chr2A_00611930_mrna',
                                             "alternative mRNA" : 'chr2A_00611930_mrna',
                                             "mismatch/match" : [2, 8, 0],
                                             "identity" : 20.0,
                                             "mismatch zones" : ([1, 5, 8, 12], [[]]),
                                             "cluster name" : "cluster 0",
                                             "reference mRNA number" : 1,
                                             "alternative mRNA number" : 1}]]
        }

    print("\n*************Testing the 'new' annotation_match function*************")
    
    for test in test_dict:
        
        print(f"\n{test} test")
        
        result = comp.annotation_match(test_dict[test][0], False, False, False)
        print(f"result : {result}\n")
        print(f"expected result : {test_dict[test][1]}\n")
        assert result == test_dict[test][1]
        
        
# test function for the 'annot_CSC.py' function 'annotation_match' (main annotation comparison function) with the 'create_strings' parameter as 'True' (uses old program version)
def test_old_annotation_match():
    
    # dictionary of inputs and expected ouputs for each test file for the 'annot_CSC.py' function 'annotation_match' (main annotation comparison function)
    test_dict = {'basic' : [cl.Cluster(name="cluster 0",
                                       loci={'ref': [locus.Locus(name='chr2A_00611930', 
                                                                    mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                                                    start=100, 
                                                                    end=299, 
                                                                    direction='direct')],
                                             'alt': [locus.Locus(name='chr2A_00611930', 
                                                                          mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                                                          start=100, 
                                                                          end=299, 
                                                                          direction='direct')]}),
                             [{"reference" : 'chr2A_00611930',
                                "reference start" : 100,
                                "reference end" : 299,
                                "alternative" : 'chr2A_00611930',
                                "alternative start" : 100,
                                "alternative end" : 299,
                                "reference mRNA" : 'chr2A_00611930_mrna',
                                "alternative mRNA" : 'chr2A_00611930_mrna',
                                "mismatch/match" : [150, 0, 0],
                                "identity" : 100.0,
                                "mismatch zones" : '?',
                                "cluster name" : "cluster 0",
                                "reference mRNA number" : 1,
                                "alternative mRNA number" : 1}]],
                 
                 'shift' : [cl.Cluster(name="cluster 0",
                                                    loci={'ref': [locus.Locus(name='chr2A_00611930', 
                                                                                 mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                                                                 start=100, 
                                                                                 end=299, 
                                                                                 direction='direct')],
                                                          'alt': [locus.Locus(name='chr2A_00611930', 
                                                                                       mRNAs={'chr2A_00611930_mrna': [100, 129, 151, 209, 240, 299]}, 
                                                                                       start=100, 
                                                                                       end=299, 
                                                                                       direction='direct')]}),
                            [{"reference" : 'chr2A_00611930',
                            "reference start" : 100,
                            "reference end" : 299,
                            "alternative" : 'chr2A_00611930',
                            "alternative start" : 100,
                            "alternative end" : 299,
                            "reference mRNA" : 'chr2A_00611930_mrna',
                            "alternative mRNA" : 'chr2A_00611930_mrna',
                            "mismatch/match" : [30, 1, 119],
                            "identity" : 20.0,
                            "mismatch zones" : '?',
                            "cluster name" : "cluster 0",
                            "reference mRNA number" : 1,
                            "alternative mRNA number" : 1}]],
                 
                 'multiple_mRNAs' : [cl.Cluster(name="cluster 0",
                                                loci={'ref': [locus.Locus(name='chr2A_00611930', 
                                                                         mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                                                         start=100, 
                                                                         end=299, 
                                                                         direction='direct')],
                                                      'alt': [locus.Locus(name='chr2A_00611930', 
                                                                         mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 179, 240, 299],
                                                                                'chr2A_00611930_mrna.2': [100, 129, 240, 299]}, 
                                                                         start=100, 
                                                                         end=299, 
                                                                         direction='direct')]}),
                                          [{"reference" : 'chr2A_00611930',
                                             "reference start" : 100,
                                             "reference end" : 299,
                                             "alternative" : 'chr2A_00611930',
                                             "alternative start" : 100,
                                             "alternative end" : 299,
                                             "reference mRNA" : 'chr2A_00611930_mrna',
                                             "alternative mRNA" : 'chr2A_00611930_mrna',
                                             "mismatch/match" : [120, 30, 0],
                                             "identity" : 80.0,
                                             "mismatch zones" : '?',
                                             "cluster name" : "cluster 0",
                                             "reference mRNA number" : 1,
                                             "alternative mRNA number" : 2}]],
                 
                 'overlapping-loci (first cluster)' : [cl.Cluster(name="cluster 0",
                                                                  loci={'ref': [locus.Locus(name='chr2A_1000', 
                                              mRNAs={'chr2A_1000_mrna': [50, 149]}, 
                                              start=50, 
                                              end=149, 
                                              direction='direct'), 
                                        locus.Locus(name='chr2A_2000', 
                                              mRNAs={'chr2A_2000_mrna': [200, 349]}, 
                                              start=200, 
                                              end=349, 
                                              direction='direct'), 
                                        locus.Locus(name='chr2A_3000', 
                                              mRNAs={'chr2A_3000_mrna': [400, 549]}, 
                                              start=400, 
                                              end=549, 
                                              direction='direct')],
                                                                        'alt': [locus.Locus(name='chr2A_1000', 
                                          mRNAs={'chr2A_1000_mrna': [100, 249]}, 
                                          start=100, 
                                          end=249, 
                                          direction='direct'), 
                                        locus.Locus(name='chr2A_2000', 
                                          mRNAs={'chr2A_2000_mrna': [300, 449]}, 
                                          start=300, 
                                          end=449, 
                                          direction='direct'), 
                                        locus.Locus(name='chr2A_3000', 
                                          mRNAs={'chr2A_3000_mrna': [500, 599]}, 
                                          start=500, 
                                          end=599, 
                                          direction='direct')]}),
                                    [{"reference" : 'chr2A_1000',
                                   "reference start" : 50,
                                   "reference end" : 149,
                                   "alternative" : 'chr2A_1000',
                                   "alternative start" : 100,
                                   "alternative end" : 249,
                                   "reference mRNA" : 'chr2A_1000_mrna',
                                   "alternative mRNA" : 'chr2A_1000_mrna',
                                   "mismatch/match" : [0, 150, 50],
                                   "identity" : 0.0,
                                   "mismatch zones" : '?',
                                   "cluster name" : "cluster 0",
                                   "reference mRNA number" : 1,
                                   "alternative mRNA number" : 1},
                                     {"reference" : 'chr2A_2000',
                                    "reference start" : 200,
                                    "reference end" : 349,
                                    "alternative" : 'chr2A_2000',
                                    "alternative start" : 300,
                                    "alternative end" : 449,
                                    "reference mRNA" : 'chr2A_2000_mrna',
                                    "alternative mRNA" : 'chr2A_2000_mrna',
                                    "mismatch/match" : [0, 200, 50],
                                    "identity" : 0.0,
                                    "mismatch zones" : '?',
                                    "cluster name" : "cluster 0",
                                    "reference mRNA number" : 1,
                                    "alternative mRNA number" : 1},
                                     {"reference" : 'chr2A_3000',
                                     "reference start" : 400,
                                     "reference end" : 549,
                                     "alternative" : 'chr2A_3000',
                                     "alternative start" : 500,
                                     "alternative end" : 599,
                                     "reference mRNA" : 'chr2A_3000_mrna',
                                     "alternative mRNA" : 'chr2A_3000_mrna',
                                     "mismatch/match" : [0, 150, 50],
                                     "identity" : 0.0,
                                     "mismatch zones" : '?',
                                     "cluster name" : "cluster 0",
                                     "reference mRNA number" : 1,
                                     "alternative mRNA number" : 1}]],
                 
                 'length_computation' : [cl.Cluster(name="cluster 0",
                                                    loci={'ref': [locus.Locus(name='chr2A_00611930', 
                                          mRNAs={'chr2A_00611930_mrna': [8, 13]}, 
                                          start=8, 
                                          end=13, 
                                          direction='direct')], 
                                                          'alt': [locus.Locus(name='chr2A_00611930', 
                                             mRNAs={'chr2A_00611930_mrna': [1, 4, 12, 13]}, 
                                             start=1, 
                                             end=13, 
                                             direction='direct')]}),
                                          [{"reference" : 'chr2A_00611930',
                                             "reference start" : 8,
                                             "reference end" : 13,
                                             "alternative" : 'chr2A_00611930',
                                             "alternative start" : 1,
                                             "alternative end" : 13,
                                             "reference mRNA" : 'chr2A_00611930_mrna',
                                             "alternative mRNA" : 'chr2A_00611930_mrna',
                                             "mismatch/match" : [2, 8, 0],
                                             "identity" : 20.0,
                                             "mismatch zones" : '?',
                                             "cluster name" : "cluster 0",
                                             "reference mRNA number" : 1,
                                             "alternative mRNA number" : 1}]]
        }

    print("\n*************Testing the 'old' annotation_match function*************")
    
    for test in test_dict:
        
        print(f"\n{test} test")
        
        result = comp.annotation_match(test_dict[test][0], True, False, False)
        print(f"result : {result}\n")
        print(f"expected result : {test_dict[test][1]}\n")
        assert result == test_dict[test][1]


# test function for the 'annot_CSC.py' function 'annotation_comparison' ('main' function of the program) with the 'create_strings' parameter as 'False' (uses new program version)
def test_new_annotation_comparison():
    
    # dictionary of inputs and expected ouputs for each test file for the 'annot_CSC.py' function 'annotation_comparison' ('main' function of the program)
    test_dict = test_dict = {
        "basic" : ["./data/tests/basic_test.gff3", 
                   "./data/tests/identical_test.gff3",
                        [[{"reference" : 'chr2A_00611930',
                           "reference start" : 100,
                           "reference end" : 299,
                           "alternative" : 'chr2A_00611930',
                           "alternative start" : 100,
                           "alternative end" : 299,
                           "reference mRNA" : 'chr2A_00611930_mrna',
                           "alternative mRNA" : 'chr2A_00611930_mrna',
                           "mismatch/match" : [150, 0, 0],
                           "identity" : 100.0,
                           "mismatch zones" : ([], [[]]),
                           "cluster name" : "cluster 0",
                           "reference mRNA number" : 1,
                           "alternative mRNA number" : 1}]]],
        
        "minus-CDS" : ["./data/tests/basic_test.gff3", 
                       "./data/tests/minus-CDS_test.gff3",
                       [[{"reference" : 'chr2A_00611930',
                          "reference start" : 100,
                          "reference end" : 299,
                          "alternative" : 'chr2A_00611930',
                          "alternative start" : 100,
                          "alternative end" : 299,
                          "reference mRNA" : 'chr2A_00611930_mrna',
                          "alternative mRNA" : 'chr2A_00611930_mrna',
                          "mismatch/match" : [90, 60, 0],
                          "identity" : 60.0,
                          "mismatch zones" : ([150, 210], [[]]),
                          "cluster name" : "cluster 0",
                          "reference mRNA number" : 1,
                          "alternative mRNA number" : 1}]]],
        
        "fusion" : ["./data/tests/basic_test.gff3", 
                    "./data/tests/fusion_test.gff3",
                    [[{"reference" : 'chr2A_00611930',
                       "reference start" : 100,
                       "reference end" : 299,
                       "alternative" : 'chr2A_00611930',
                       "alternative start" : 100,
                       "alternative end" : 299,
                       "reference mRNA" : 'chr2A_00611930_mrna',
                       "alternative mRNA" : 'chr2A_00611930_mrna',
                       "mismatch/match" : [30, 20, 120],
                       "identity" : 17.6,
                       "mismatch zones" : ([130, 150], [[150, 209], [240, 299]]),
                       "cluster name" : "cluster 0",
                       "reference mRNA number" : 1,
                       "alternative mRNA number" : 1}]]],
        
        "shift" : ["./data/tests/basic_test.gff3", 
                   "./data/tests/shift_test.gff3",
                   [[{"reference" : 'chr2A_00611930',
                      "reference start" : 100,
                      "reference end" : 299,
                      "alternative" : 'chr2A_00611930',
                      "alternative start" : 100,
                      "alternative end" : 299,
                      "reference mRNA" : 'chr2A_00611930_mrna',
                      "alternative mRNA" : 'chr2A_00611930_mrna',
                      "mismatch/match" : [30, 1, 119],
                      "identity" : 20.0,
                      "mismatch zones" : ([150, 151], [[151, 209], [240, 299]]),
                      "cluster name" : "cluster 0",
                      "reference mRNA number" : 1,
                      "alternative mRNA number" : 1}]]],
        
        "reverse-basic" : ["./data/tests/basic_test.gff3", 
                     "./data/tests/reverse_test.gff3",
                     [[{"reference" : '~',
                        "reference start" : '_',
                        "reference end" : '_',
                        "alternative" : 'chr2A_00611930',
                        "alternative start" : 100,
                        "alternative end" : 299,
                        "reference mRNA" : '_',
                        "alternative mRNA" : '_',
                        "mismatch/match" : [],
                        "identity" : 0.0,
                        "mismatch zones" : '_',
                        "cluster name" : "cluster 0",
                        "reference mRNA number" : '_',
                        "alternative mRNA number" : 1}],
                      [{"reference" : 'chr2A_00611930',
                        "reference start" : 100,
                        "reference end" : 299,
                        "alternative" : '~',
                        "alternative start" : '_',
                        "alternative end" : '_',
                        "reference mRNA" : '_',
                        "alternative mRNA" : '_',
                        "mismatch/match" : [],
                        "identity" : 0.0,
                        "mismatch zones" : '_',
                        "cluster name" : "cluster 1",
                        "reference mRNA number" : 1,
                        "alternative mRNA number" : '_'}]]],
        
        "reverse-reverse" : ["./data/tests/reverse_test.gff3", 
                     "./data/tests/reverse_test.gff3",
                     [[{"reference" : 'chr2A_00611930',
                        "reference start" : 100,
                        "reference end" : 299,
                        "alternative" : 'chr2A_00611930',
                        "alternative start" : 100,
                        "alternative end" : 299,
                        "reference mRNA" : 'chr2A_00611930_mrna',
                        "alternative mRNA" : 'chr2A_00611930_mrna',
                        "mismatch/match" : [150, 0, 0],
                        "identity" : 100.0,
                        "mismatch zones" : ([], []),
                        "cluster name" : "cluster 0",
                        "reference mRNA number" : 1,
                        "alternative mRNA number" : 1}]]],     
        
        "diff-start-before" : ["./data/tests/basic_test.gff3",
                               "./data/tests/diff-start-before_test.gff3",
                               [[{"reference" : 'chr2A_00611930',
                                  "reference start" : 100,
                                  "reference end" : 299,
                                  "alternative" : 'chr2A_00611930',
                                  "alternative start" : 40,
                                  "alternative end" : 239,
                                  "reference mRNA" : 'chr2A_00611930_mrna',
                                  "alternative mRNA" : 'chr2A_00611930_mrna',
                                  "mismatch/match" : [30, 180, 30],
                                  "identity" : 12.5,
                                  "mismatch zones" : ([40, 70, 90, 100, 130, 180, 210, 300], [[100, 129]]),
                                  "cluster name" : "cluster 0",
                                  "reference mRNA number" : 1,
                                  "alternative mRNA number" : 1}]]],
        
        "diff-start-after" : ["./data/tests/basic_test.gff3",
                              "./data/tests/diff-start-after_test.gff3",
                              [[{"reference" : 'chr2A_00611930',
                                 "reference start" : 100,
                                 "reference end" : 299,
                                 "alternative" : 'chr2A_00611930',
                                 "alternative start" : 160,
                                 "alternative end" : 359,
                                 "reference mRNA" : 'chr2A_00611930_mrna',
                                 "alternative mRNA" : 'chr2A_00611930_mrna',
                                 "mismatch/match" : [30, 180, 30],
                                 "identity" : 12.5,
                                 "mismatch zones" : ([100, 130, 150, 160, 190, 240, 270, 360], [[160, 189]]),
                                 "cluster name" : "cluster 0",
                                 "reference mRNA number" : 1,
                                 "alternative mRNA number" : 1}]]],
        
        "basic-2-loci" : ["./data/tests/basic-2-loci_test.gff3",
                          "./data/tests/identical-2-loci_test.gff3",
                              [[{"reference" : 'chr2A_00611930',
                                 "reference start" : 100,
                                 "reference end" : 299,
                                 "alternative" : 'chr2A_00611930',
                                 "alternative start" : 100,
                                 "alternative end" : 299,
                                 "reference mRNA" : 'chr2A_00611930_mrna',
                                 "alternative mRNA" : 'chr2A_00611930_mrna',
                                 "mismatch/match" : [150, 0, 0],
                                 "identity" : 100.0,
                                 "mismatch zones" : ([], [[]]),
                                 "cluster name" : "cluster 0",
                                 "reference mRNA number" : 1,
                                 "alternative mRNA number" : 1}],
                               [{"reference" : 'chr2A_00620000',
                                  "reference start" : 600,
                                  "reference end" : 899,
                                  "alternative" : 'chr2A_00620000',
                                  "alternative start" : 600,
                                  "alternative end" : 899,
                                  "reference mRNA" : 'chr2A_00620000_mrna',
                                  "alternative mRNA" : 'chr2A_00620000_mrna',
                                  "mismatch/match" : [200, 0, 0],
                                  "identity" : 100.0,
                                  "mismatch zones" : ([], [[]]),
                                  "cluster name" : "cluster 2",
                                  "reference mRNA number" : 1,
                                  "alternative mRNA number" : 1}]]],
        
        "minus-loci" : ["./data/tests/basic_test.gff3",
                          "./data/tests/identical-2-loci_test.gff3",
                              [[{"reference" : 'chr2A_00611930',
                                 "reference start" : 100,
                                 "reference end" : 299,
                                 "alternative" : 'chr2A_00611930',
                                 "alternative start" : 100,
                                 "alternative end" : 299,
                                 "reference mRNA" : 'chr2A_00611930_mrna',
                                 "alternative mRNA" : 'chr2A_00611930_mrna',
                                 "mismatch/match" : [150, 0, 0],
                                 "identity" : 100.0,
                                 "mismatch zones" : ([], [[]]),
                                 "cluster name" : "cluster 0",
                                 "reference mRNA number" : 1,
                                 "alternative mRNA number" : 1}],
                               [{"reference" : '~',
                                  "reference start" : '_',
                                  "reference end" : '_',
                                  "alternative" : 'chr2A_00620000',
                                  "alternative start" : 600,
                                  "alternative end" : 899,
                                  "reference mRNA" : '_',
                                  "alternative mRNA" : '_',
                                  "mismatch/match" : [],
                                  "identity" : 0.0,
                                  "mismatch zones" : '_',
                                  "cluster name" : "cluster 2",
                                  "reference mRNA number" : '_',
                                  "alternative mRNA number" : 1}]]],
        
        "overlapping-loci" : ["./data/tests/overlapping-loci_test.gff3",
                              "./data/tests/overlapping-loci-alt_test.gff3",
                              [[{'reference': 'chr2A_1000',
                               'reference start': 50,
                               'reference end': 149,
                               'alternative': 'chr2A_1000',
                               'alternative start': 100,
                               'alternative end': 249,
                               "reference mRNA" : 'chr2A_1000_mrna',
                               "alternative mRNA" : 'chr2A_1000_mrna',
                               'mismatch/match': [0, 150, 50],
                               'identity': 0.0,
                               'mismatch zones': ([50, 100, 150, 250], [[100, 149]]),
                               "cluster name" : "cluster 0",
                               "reference mRNA number" : 1,
                               "alternative mRNA number" : 1},
                                {'reference': 'chr2A_2000',
                                 'reference start': 200,
                                 'reference end': 349,
                                 'alternative': 'chr2A_2000',
                                 'alternative start': 300,
                                 'alternative end': 449,
                                 "reference mRNA" : 'chr2A_2000_mrna',
                                 "alternative mRNA" : 'chr2A_2000_mrna',
                                 'mismatch/match': [0, 200, 50],
                                 'identity': 0.0,
                                 'mismatch zones': ([200, 300, 350, 450], [[300, 349]]),
                                 "cluster name" : "cluster 0",
                                 "reference mRNA number" : 1,
                                 "alternative mRNA number" : 1},
                                {'reference': 'chr2A_3000',
                                 'reference start': 400, 
                                 'reference end': 549,
                                 'alternative': 'chr2A_3000', 
                                 'alternative start': 500,
                                 'alternative end': 599, 
                                 "reference mRNA" : 'chr2A_3000_mrna',
                                 "alternative mRNA" : 'chr2A_3000_mrna',
                                 'mismatch/match': [0, 150, 50], 
                                 'identity': 0.0, 
                                 'mismatch zones': ([400, 500, 550, 600], [[500, 549]]),
                                 "cluster name" : "cluster 0",
                                 "reference mRNA number" : 1,
                                 "alternative mRNA number" : 1}],
                               [{'reference': 'chr2A_4000', 
                                 'reference start': 650,
                                 'reference end': 699, 
                                 'alternative': 'chr2A_4000',
                                 'alternative start': 650, 
                                 'alternative end': 699, 
                                 "reference mRNA" : 'chr2A_4000_mrna',
                                 "alternative mRNA" : 'chr2A_4000_mrna',
                                 'mismatch/match': [50, 0, 0], 
                                 'identity': 100.0, 
                                 'mismatch zones': ([], [[]]),
                                 "cluster name" : "cluster 6",
                                 "reference mRNA number" : 1,
                                 "alternative mRNA number" : 1}], 
                               [{'reference': 'chr2A_5000',
                                 'reference start': 750, 
                                 'reference end': 799,
                                 'alternative': 'chr2A_5000', 
                                 'alternative start': 750,
                                 'alternative end': 779, 
                                 "reference mRNA" : 'chr2A_5000_mrna',
                                 "alternative mRNA" : 'chr2A_5000_mrna',
                                 'mismatch/match': [30, 20, 0], 
                                 'identity': 60.0,
                                 'mismatch zones': ([780, 800], [[]]),
                                 "cluster name" : "cluster 8",
                                 "reference mRNA number" : 1,
                                 "alternative mRNA number" : 1}, 
                                {'reference': '~',
                                 'reference start': '_', 
                                 'reference end': '_',
                                 'alternative': 'chr2A_6000',
                                 'alternative start': 790, 
                                 'alternative end': 849,
                                 "reference mRNA" : '_',
                                 "alternative mRNA" : 'chr2A_6000_mrna',
                                 'mismatch/match': [], 
                                 'identity': '_', 
                                 'mismatch zones': '_',
                                 "cluster name" : "cluster 8",
                                 "reference mRNA number" : '_',
                                 "alternative mRNA number" : 1}]]],
        
        'length_computation' : ["./data/tests/length_computation_ref_test.gff3", 
                                "./data/tests/length_computation_alt_test.gff3",
                                 [[{"reference" : 'chr2A_00611930',
                                    "reference start" : 8,
                                    "reference end" : 13,
                                    "alternative" : 'chr2A_00611930',
                                    "alternative start" : 1,
                                    "alternative end" : 13,
                                    "reference mRNA" : 'chr2A_00611930_mrna',
                                    "alternative mRNA" : 'chr2A_00611930_mrna',
                                    "mismatch/match" : [2, 8, 0],
                                    "identity" : 20.0,
                                    "mismatch zones" : ([1, 5, 8, 12], [[]]),
                                    "cluster name" : "cluster 0",
                                    "reference mRNA number" : 1,
                                    "alternative mRNA number" : 1}]]]
        
        }

    print("\n*************Testing the 'new' annotation_comparison function*************")
    
    for test in test_dict:
        
        print(f"\n{test} test")
        
        result = annot_CSC.annotation_comparison(test_dict[test][0], test_dict[test][1], False, False, False)
        print(f"result : {result}\n")
        print(f"expected result : {test_dict[test][2]}\n")
        assert result == test_dict[test][2]
        
        
# test function for the 'annot_CSC.py' function 'annotation_comparison' ('main' function of the program) with the 'create_strings' parameter as 'True' (uses old program version)
def test_old_annotation_comparison():
    
    # dictionary of inputs and expected ouputs for each test file for the 'annot_CSC.py' function 'annotation_comparison' ('main' function of the program)
    test_dict = test_dict = {
        "basic" : ["./data/tests/basic_test.gff3", 
                   "./data/tests/identical_test.gff3",
                        [[{"reference" : 'chr2A_00611930',
                           "reference start" : 100,
                           "reference end" : 299,
                           "alternative" : 'chr2A_00611930',
                           "alternative start" : 100,
                           "alternative end" : 299,
                           "reference mRNA" : 'chr2A_00611930_mrna',
                           "alternative mRNA" : 'chr2A_00611930_mrna',
                           "mismatch/match" : [150, 0, 0],
                           "identity" : 100.0,
                           "mismatch zones" : '?',
                           "cluster name" : "cluster 0",
                           "reference mRNA number" : 1,
                           "alternative mRNA number" : 1}]]],
        
        "minus-CDS" : ["./data/tests/basic_test.gff3", 
                       "./data/tests/minus-CDS_test.gff3",
                       [[{"reference" : 'chr2A_00611930',
                          "reference start" : 100,
                          "reference end" : 299,
                          "alternative" : 'chr2A_00611930',
                          "alternative start" : 100,
                          "alternative end" : 299,
                          "reference mRNA" : 'chr2A_00611930_mrna',
                          "alternative mRNA" : 'chr2A_00611930_mrna',
                          "mismatch/match" : [90, 60, 0],
                          "identity" : 60.0,
                          "mismatch zones" : '?',
                          "cluster name" : "cluster 0",
                          "reference mRNA number" : 1,
                          "alternative mRNA number" : 1}]]],
        
        "fusion" : ["./data/tests/basic_test.gff3", 
                    "./data/tests/fusion_test.gff3",
                    [[{"reference" : 'chr2A_00611930',
                       "reference start" : 100,
                       "reference end" : 299,
                       "alternative" : 'chr2A_00611930',
                       "alternative start" : 100,
                       "alternative end" : 299,
                       "reference mRNA" : 'chr2A_00611930_mrna',
                       "alternative mRNA" : 'chr2A_00611930_mrna',
                       "mismatch/match" : [30, 20, 120],
                       "identity" : 17.6,
                       "mismatch zones" : '?',
                       "cluster name" : "cluster 0",
                       "reference mRNA number" : 1,
                       "alternative mRNA number" : 1}]]],
        
        "shift" : ["./data/tests/basic_test.gff3", 
                   "./data/tests/shift_test.gff3",
                   [[{"reference" : 'chr2A_00611930',
                      "reference start" : 100,
                      "reference end" : 299,
                      "alternative" : 'chr2A_00611930',
                      "alternative start" : 100,
                      "alternative end" : 299,
                      "reference mRNA" : 'chr2A_00611930_mrna',
                      "alternative mRNA" : 'chr2A_00611930_mrna',
                      "mismatch/match" : [30, 1, 119],
                      "identity" : 20.0,
                      "mismatch zones" : '?',
                      "cluster name" : "cluster 0",
                      "reference mRNA number" : 1,
                      "alternative mRNA number" : 1}]]],
        
        "reverse-basic" : ["./data/tests/basic_test.gff3", 
                     "./data/tests/reverse_test.gff3",
                     [[{"reference" : '~',
                        "reference start" : '_',
                        "reference end" : '_',
                        "alternative" : 'chr2A_00611930',
                        "alternative start" : 100,
                        "alternative end" : 299,
                        "reference mRNA" : '_',
                        "alternative mRNA" : '_',
                        "mismatch/match" : [],
                        "identity" : 0.0,
                        "mismatch zones" : '_',
                        "cluster name" : "cluster 0",
                        "reference mRNA number" : '_',
                        "alternative mRNA number" : 1}],
                     [{"reference" : 'chr2A_00611930',
                        "reference start" : 100,
                        "reference end" : 299,
                        "alternative" : '~',
                        "alternative start" : '_',
                        "alternative end" : '_',
                        "reference mRNA" : '_',
                        "alternative mRNA" : '_',
                        "mismatch/match" : [],
                        "identity" : 0.0,
                        "mismatch zones" : '_',
                        "cluster name" : "cluster 1",
                        "reference mRNA number" : 1,
                        "alternative mRNA number" : '_'}]]],
        
        "reverse-reverse" : ["./data/tests/reverse_test.gff3", 
                     "./data/tests/reverse_test.gff3",
                     [[{"reference" : 'chr2A_00611930',
                        "reference start" : 100,
                        "reference end" : 299,
                        "alternative" : 'chr2A_00611930',
                        "alternative start" : 100,
                        "alternative end" : 299,
                        "reference mRNA" : 'chr2A_00611930_mrna',
                        "alternative mRNA" : 'chr2A_00611930_mrna',
                        "mismatch/match" : [150, 0, 0],
                        "identity" : 100.0,
                        "mismatch zones" : '?',
                        "cluster name" : "cluster 0",
                        "reference mRNA number" : 1,
                        "alternative mRNA number" : 1}]]],        
        
        "diff-start-before" : ["./data/tests/basic_test.gff3",
                               "./data/tests/diff-start-before_test.gff3",
                               [[{"reference" : 'chr2A_00611930',
                                  "reference start" : 100,
                                  "reference end" : 299,
                                  "alternative" : 'chr2A_00611930',
                                  "alternative start" : 40,
                                  "alternative end" : 239,
                                  "reference mRNA" : 'chr2A_00611930_mrna',
                                  "alternative mRNA" : 'chr2A_00611930_mrna',
                                  "mismatch/match" : [30, 180, 30],
                                  "identity" : 12.5,
                                  "mismatch zones" : '?',
                                  "cluster name" : "cluster 0",
                                  "reference mRNA number" : 1,
                                  "alternative mRNA number" : 1}]]],
        
        "diff-start-after" : ["./data/tests/basic_test.gff3",
                              "./data/tests/diff-start-after_test.gff3",
                              [[{"reference" : 'chr2A_00611930',
                                 "reference start" : 100,
                                 "reference end" : 299,
                                 "alternative" : 'chr2A_00611930',
                                 "alternative start" : 160,
                                 "alternative end" : 359,
                                 "reference mRNA" : 'chr2A_00611930_mrna',
                                 "alternative mRNA" : 'chr2A_00611930_mrna',
                                 "mismatch/match" : [30, 180, 30],
                                 "identity" : 12.5,
                                 "mismatch zones" : '?',
                                 "cluster name" : "cluster 0",
                                 "reference mRNA number" : 1,
                                 "alternative mRNA number" : 1}]]],
        
        "basic-2-loci" : ["./data/tests/basic-2-loci_test.gff3",
                          "./data/tests/identical-2-loci_test.gff3",
                              [[{"reference" : 'chr2A_00611930',
                                 "reference start" : 100,
                                 "reference end" : 299,
                                 "alternative" : 'chr2A_00611930',
                                 "alternative start" : 100,
                                 "alternative end" : 299,
                                 "reference mRNA" : 'chr2A_00611930_mrna',
                                 "alternative mRNA" : 'chr2A_00611930_mrna',
                                 "mismatch/match" : [150, 0, 0],
                                 "identity" : 100.0,
                                 "mismatch zones" : '?',
                                 "cluster name" : "cluster 0",
                                 "reference mRNA number" : 1,
                                 "alternative mRNA number" : 1}],
                               [{"reference" : 'chr2A_00620000',
                                  "reference start" : 600,
                                  "reference end" : 899,
                                  "alternative" : 'chr2A_00620000',
                                  "alternative start" : 600,
                                  "alternative end" : 899,
                                  "reference mRNA" : 'chr2A_00620000_mrna',
                                  "alternative mRNA" : 'chr2A_00620000_mrna',
                                  "mismatch/match" : [200, 0, 0],
                                  "identity" : 100.0,
                                  "mismatch zones" : '?',
                                  "cluster name" : "cluster 2",
                                  "reference mRNA number" : 1,
                                  "alternative mRNA number" : 1}]]],
        
        "minus-loci" : ["./data/tests/basic_test.gff3",
                          "./data/tests/identical-2-loci_test.gff3",
                              [[{"reference" : 'chr2A_00611930',
                                 "reference start" : 100,
                                 "reference end" : 299,
                                 "alternative" : 'chr2A_00611930',
                                 "alternative start" : 100,
                                 "alternative end" : 299,
                                 "reference mRNA" : 'chr2A_00611930_mrna',
                                 "alternative mRNA" : 'chr2A_00611930_mrna',
                                 "mismatch/match" : [150, 0, 0],
                                 "identity" : 100.0,
                                 "mismatch zones" : '?',
                                 "cluster name" : "cluster 0",
                                 "reference mRNA number" : 1,
                                 "alternative mRNA number" : 1}],
                               [{"reference" : '~',
                                  "reference start" : '_',
                                  "reference end" : '_',
                                  "alternative" : 'chr2A_00620000',
                                  "alternative start" : 600,
                                  "alternative end" : 899,
                                  "reference mRNA" : '_',
                                  "alternative mRNA" : '_',
                                  "mismatch/match" : [],
                                  "identity" : 0.0,
                                  "mismatch zones" : '_',
                                  "cluster name" : "cluster 2",
                                  "reference mRNA number" : '_',
                                  "alternative mRNA number" : 1}]]],
        
        "overlapping-loci" : ["./data/tests/overlapping-loci_test.gff3",
                              "./data/tests/overlapping-loci-alt_test.gff3",
                              [[{'reference': 'chr2A_1000',
                               'reference start': 50,
                               'reference end': 149,
                               'alternative': 'chr2A_1000',
                               'alternative start': 100,
                               'alternative end': 249,
                               "reference mRNA" : 'chr2A_1000_mrna',
                               "alternative mRNA" : 'chr2A_1000_mrna',
                               'mismatch/match': [0, 150, 50],
                               'identity': 0.0,
                               'mismatch zones': '?',
                               "cluster name" : "cluster 0",
                               "reference mRNA number" : 1,
                               "alternative mRNA number" : 1},
                                {'reference': 'chr2A_2000',
                                 'reference start': 200,
                                 'reference end': 349,
                                 'alternative': 'chr2A_2000',
                                 'alternative start': 300,
                                 'alternative end': 449,
                                 "reference mRNA" : 'chr2A_2000_mrna',
                                 "alternative mRNA" : 'chr2A_2000_mrna',
                                 'mismatch/match': [0, 200, 50],
                                 'identity': 0.0,
                                 'mismatch zones': '?',
                                 "cluster name" : "cluster 0",
                                 "reference mRNA number" : 1,
                                 "alternative mRNA number" : 1},
                                {'reference': 'chr2A_3000',
                                 'reference start': 400, 
                                 'reference end': 549,
                                 'alternative': 'chr2A_3000', 
                                 'alternative start': 500,
                                 'alternative end': 599, 
                                 "reference mRNA" : 'chr2A_3000_mrna',
                                 "alternative mRNA" : 'chr2A_3000_mrna',
                                 'mismatch/match': [0, 150, 50], 
                                 'identity': 0.0, 
                                 'mismatch zones': '?',
                                 "cluster name" : "cluster 0",
                                 "reference mRNA number" : 1,
                                 "alternative mRNA number" : 1}],
                               [{'reference': 'chr2A_4000', 
                                 'reference start': 650,
                                 'reference end': 699, 
                                 'alternative': 'chr2A_4000',
                                 'alternative start': 650, 
                                 'alternative end': 699, 
                                 "reference mRNA" : 'chr2A_4000_mrna',
                                 "alternative mRNA" : 'chr2A_4000_mrna',
                                 'mismatch/match': [50, 0, 0], 
                                 'identity': 100.0, 
                                 'mismatch zones': '?',
                                 "cluster name" : "cluster 6",
                                 "reference mRNA number" : 1,
                                 "alternative mRNA number" : 1}], 
                               [{'reference': 'chr2A_5000',
                                 'reference start': 750, 
                                 'reference end': 799,
                                 'alternative': 'chr2A_5000', 
                                 'alternative start': 750,
                                 'alternative end': 779, 
                                 "reference mRNA" : 'chr2A_5000_mrna',
                                 "alternative mRNA" : 'chr2A_5000_mrna',
                                 'mismatch/match': [30, 20, 0], 
                                 'identity': 60.0,
                                 'mismatch zones': '?',
                                 "cluster name" : "cluster 8",
                                 "reference mRNA number" : 1,
                                 "alternative mRNA number" : 1}, 
                                {'reference': '~',
                                 'reference start': '_', 
                                 'reference end': '_',
                                 'alternative': 'chr2A_6000',
                                 'alternative start': 790, 
                                 'alternative end': 849,
                                 "reference mRNA" : '_',
                                 "alternative mRNA" : 'chr2A_6000_mrna',
                                 'mismatch/match': [], 
                                 'identity': '_', 
                                 'mismatch zones': '_',
                                 "cluster name" : "cluster 8",
                                 "reference mRNA number" : '_',
                                 "alternative mRNA number" : 1}]]],
        
        'length_computation' : ["./data/tests/length_computation_ref_test.gff3", 
                                "./data/tests/length_computation_alt_test.gff3",
                                 [[{"reference" : 'chr2A_00611930',
                                    "reference start" : 8,
                                    "reference end" : 13,
                                    "alternative" : 'chr2A_00611930',
                                    "alternative start" : 1,
                                    "alternative end" : 13,
                                    "reference mRNA" : 'chr2A_00611930_mrna',
                                    "alternative mRNA" : 'chr2A_00611930_mrna',
                                    "mismatch/match" : [2, 8, 0],
                                    "identity" : 20.0,
                                    "mismatch zones" : '?',
                                    "cluster name" : "cluster 0",
                                    "reference mRNA number" : 1,
                                    "alternative mRNA number" : 1}]]]
        
        }

    print("\n*************Testing the 'old' annotation_comparison function*************")
    
    for test in test_dict:
        
        print(f"\n{test} test")
        
        result = annot_CSC.annotation_comparison(test_dict[test][0], test_dict[test][1], False, False, True)
        print(f"result : {result}\n")
        print(f"expected result : {test_dict[test][2]}\n")
        assert result == test_dict[test][2]
        
        
# test function for the 'comparison.py' function 'compute_matches_mismatches_EI_RF' (CDS intervals comparison function)
def test_compute_matches_mismatches_EI_RF():
    
    # dictionary of inputs and expected ouputs for the function
    test_dict = {
        "identical" : [[100, 129, 150, 209, 240, 299],
                       iu.OrderedIntervals(intervals=[100, 129, 150, 209, 240, 299], include_ub=True), 
                       [100, 129, 150, 209, 240, 299],
                       (150, 0, 0, [], [[]])],
        
        "fusion" : [[100, 129, 150, 209, 240, 299],
                       iu.OrderedIntervals(intervals=[100, 209, 240, 299], include_ub=True), 
                       [100, 209, 240, 299],
                       (30, 20, 120, [130, 150], [[150, 209], [240, 299]])],
        
        "shift" : [[100, 129, 150, 209, 240, 299],
                       iu.OrderedIntervals(intervals=[100, 129, 150, 209, 240, 299], include_ub=True), 
                       [100, 129, 151, 209, 240, 299],
                       (30, 1, 119, [150, 151], [[151, 209], [240, 299]])],
        
        "reverse-reverse" : [[0, 59, 90, 149, 170, 199],
                             iu.OrderedIntervals(intervals=[0, 60, 90, 150, 170, 200], include_ub=True), 
                             [0, 59, 90, 149, 170, 199],
                             (150, 0, 0, [], [[]])],
        
        "diff-start-before" : [[100, 129, 150, 209, 240, 299],
                       iu.OrderedIntervals(intervals=[100, 130, 150, 210, 240, 300], include_ub=True), 
                       [40, 69, 90, 149, 180, 239],
                       (30, 180, 30, [40, 70, 90, 100, 130, 180, 210, 300], [[100, 129]])],
        
        "diff-start-after" : [[100, 129, 150, 209, 240, 299],
                       iu.OrderedIntervals(intervals=[100, 130, 150, 210, 240, 300], include_ub=True), 
                       [160, 189, 210, 269, 300, 359],
                       (30, 180, 30, [100, 130, 150, 160, 190, 240, 270, 360], [[160, 189]])],
        
        'overlapping-loci (first cluster/locus)' : [[50, 149],
                       iu.OrderedIntervals(intervals=[50, 150], include_ub=True), 
                       [100, 249],
                       (0, 150, 50, [50, 100, 150, 250], [[100, 149]])],
        
        'overlapping-loci (first cluster / 2nd locus)' : [[200, 349],
                       iu.OrderedIntervals(intervals=[200, 350], include_ub=True), 
                       [100, 249],
                       (0, 200, 50, [100, 200, 250, 350], [[200, 249]])],
                         
        "length_computation" : [[8, 13],
                       iu.OrderedIntervals(intervals=[8, 14], include_ub=True), 
                       [1, 4, 12, 13],
                       (2, 8, 0, [1, 5, 8, 12], [[]])]
        }
    
    print("\n*************Testing the compute_matches_mismatches_EI_RF function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        result = comp.compute_matches_mismatches_EI_RF(test_dict[test][0], test_dict[test][1], test_dict[test][2], False, False)
        print(f"result : {result}\n")
        print(f"expected result : {test_dict[test][3]}\n")
        assert result == test_dict[test][3]



