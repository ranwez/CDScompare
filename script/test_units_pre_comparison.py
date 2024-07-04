#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 11:00:43 2024

@author: vetea
"""


import os, sys

# code adapted from https://csatlas.com/python-import-file-module/#import_a_file_in_a_different_directory
script_dir = os.path.dirname( __file__ )
sys.path.append( script_dir )
import pre_comparison as pc
import locus


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
        
        "reverse" : [[0, 59, 90, 149, 170, 199],
                     [0, 59, 90, 149, 170, 199],
                     [1, 1, 1]],
        
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
        

# test function for the 'pre_comparison.py' function 'annotation_sort' (locus list creation and sorting function)
def test_annotation_sort():
    
    # dictionary of inputs and expected ouputs for each test file for the 'pre_comparison.py' function 'annotation_sort' (locus list creation and sorting function)
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
        
        result = pc.annotation_sort(test_dict[test][0], test_dict[test][1], False)
        print(f"result : {result}\n")
        print(f"expected result : {test_dict[test][2]}\n")
        assert result == test_dict[test][2]
        
        
# test function for the 'pre_comparison.py' function 'construct_clusters' (overlapping loci grouping function)
def test_construct_clusters():
    
    # dictionary of inputs and expected ouputs for each test file for the 'pre_comparison.py' function 'construct_clusters' (overlapping loci grouping function)
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
        
        result = pc.construct_clusters(test_dict[test][0], test_dict[test][1], test_dict[test][2], False)
        for cluster_id, cluster in test_dict[test][3].items():
            result_details = result[cluster_id].get_details()
            print(f"expected mRNAs: {test_dict[test][3][cluster_id]}")
            print(f"result mRNAs: {result_details}")
            assert result_details == cluster
            
            
# test function for the 'pre_comparison.py' function 'create_vectors' (structure string creation function)
def test_create_vectors():
    
    # dictionary of inputs and expected ouputs for each test file for the 'pre_comparison.py' function 'create_vectors' (structure string creation function)
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
        
        result = pc.create_vectors(test_dict[test][0], False)
        print(f"result : {result}\n")
        print(test_dict[test][1])
        assert result == test_dict[test][1]