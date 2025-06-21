#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 10:53:48 2024

@author: vetea
"""


import os, sys

script_dir = os.path.dirname( __file__ )
script_dir = "/".join(script_dir.split("/")[:-1]) + "/python_util/"
sys.path.append( script_dir )

import comparison as comp
import locus
import cluster as cl
import intervals_utils as iu
from comparison import MrnaMatchInfo, MismatchInfo

# test function for the 'comparison.py' class function 'reverse_coord' (comparison mismatch zones inversion function)
def test_reverse_coord():
    
    # dictionary of inputs and expected ouputs for the method
    test_dict = {
        
        "reverse-reverse" : [([], []),
                             299,
                             ([], [])],
        
        "reverse-reverse-modified" : [([80, 90], [90, 149, 170, 199]),
                                      299,
                                      ([209, 219], [100, 129, 150, 209])]
        
        }
    
    print("\n*************Testing the reverse_coord function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        result = comp.reverse_coord(test_dict[test][0], test_dict[test][1])
        print(f"result : {result}\n")
        print(test_dict[test][2])
        assert result == test_dict[test][2]
        
        
# test function for the 'comparison.py' function 'compare_loci' (new locus comparison function)
def test_compare_loci():
    
    # dictionary of inputs and expected ouputs for each test file for the 'comparison.py' function 'compare_loci' (new locus comparison function)
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
                       ([150, 0, 0], 100.0, ([], []), 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
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
                       ([90, 60, 0], 60.0, ([150, 210], []), 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
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
                    ([30, 20, 120], 17.6, ([130, 150], [150, 209, 240, 299]), 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
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
                   ([30, 1, 119], 20.0, ([150, 151], [151, 209, 240, 299]), 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
        "reverse" : [locus.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [0, 59, 90, 149, 170, 199]}, 
                                 start=100, 
                                 end=299, 
                                 direction='reverse'),
                     locus.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [0, 59, 90, 149, 170, 199]}, 
                                 start=100, 
                                 end=299, 
                                 direction='reverse'),
                     ([150, 0, 0], 100.0, ([], []), 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
        "reverse-modified" : [locus.Locus(name='chr2A_00611930', 
                                          mRNAs={'chr2A_00611930_mrna': [0, 59, 90, 149, 170, 199]}, 
                                          start=100, 
                                          end=299, 
                                          direction='reverse'),
                              locus.Locus(name='chr2A_00611930', 
                                          mRNAs={'chr2A_00611930_mrna': [0, 59, 80, 149, 170, 199]}, 
                                          start=100, 
                                          end=299, 
                                          direction='reverse'),
                              ([60, 10, 90], 37.5, ([80, 90], [90, 149, 170, 199]), 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
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
                               ([30, 180, 30], 12.5, ([40, 70, 90, 100, 130, 180, 210, 300], [100, 129]), 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
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
                              ([30, 180, 30], 12.5, ([100, 130, 150, 160, 190, 240, 270, 360], [160, 189]), 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
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
                       ([120, 30, 0], 80.0, ([180, 210], []), 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
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
                                         ([200, 0, 0], 100.0, ([], []), 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')],
        
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
                       ([2, 8, 0], 20.0, ([1, 5, 8, 12], []), 'chr2A_00611930_mrna', 'chr2A_00611930_mrna')]
        
        }
    
    print("\n*************Testing the compare_loci function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        result = comp.compare_loci(test_dict[test][0], test_dict[test][1], False)
        print(f"Expected result = {test_dict[test][2]}")
        print(f"Result = {result}")
        assert result == test_dict[test][2]
        
           
        
# test function for the 'comparison.py' function 'annotation_match' (main annotation comparison function) with the 'create_strings' parameter as 'False' (uses new program version)
def test_new_annotation_match():
    
    # dictionary of inputs and expected ouputs for each test file for the 'comparison.py' function 'annotation_match' (main annotation comparison function)
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
                                "mismatch zones" : ([], []),
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
                            "mismatch zones" : ([150, 151], [151, 209, 240, 299]),
                            "cluster name" : "cluster 0",
                            "reference mRNA number" : 1,
                            "alternative mRNA number" : 1}]],
                 
                 'reverse' : [cl.Cluster(name="cluster 0",
                                         end=299,
                                                    loci={'ref': [locus.Locus(name='chr2A_00611930', 
                                                                             mRNAs={'chr2A_00611930_mrna': [0, 59, 90, 149, 170, 199]}, 
                                                                             start=100, 
                                                                             end=299, 
                                                                             direction='reverse')],
                                                          'alt': [locus.Locus(name='chr2A_00611930', 
                                                                                   mRNAs={'chr2A_00611930_mrna': [0, 59, 90, 149, 170, 199]}, 
                                                                                   start=100, 
                                                                                   end=299, 
                                                                                   direction='reverse')]}),
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
                                             "mismatch zones" : ([], []),
                                             "cluster name" : "cluster 0",
                                             "reference mRNA number" : 1,
                                             "alternative mRNA number" : 1}]],
                 
                 'reverse-modified' : [cl.Cluster(name="cluster 0",
                                                  end=299,
                                                    loci={'ref': [locus.Locus(name='chr2A_00611930', 
                                                                             mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                                                             start=100, 
                                                                             end=299, 
                                                                             direction='reverse')],
                                                          'alt': [locus.Locus(name='chr2A_00611930', 
                                                                                   mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 219, 240, 299]}, 
                                                                                   start=100, 
                                                                                   end=299, 
                                                                                   direction='reverse')]}),
                                          [{"reference" : 'chr2A_00611930',
                                             "reference start" : 100,
                                             "reference end" : 299,
                                             "alternative" : 'chr2A_00611930',
                                             "alternative start" : 100,
                                             "alternative end" : 299,
                                             "reference mRNA" : 'chr2A_00611930_mrna',
                                             "alternative mRNA" : 'chr2A_00611930_mrna',
                                             "mismatch/match" : [60, 10, 90],
                                             "identity" : 37.5,
                                             "mismatch zones" : ([209, 219], [100, 129, 150, 209]),
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
                                             "mismatch zones" : ([180, 210], []),
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
                                   "mismatch zones" : ([50, 100, 150, 250], [100, 149]),
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
                                    "mismatch zones" : ([200, 300, 350, 450], [300, 349]),
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
                                     "mismatch zones" : ([400, 500, 550, 600], [500, 549]),
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
                                             "mismatch zones" : ([1, 5, 8, 12], []),
                                             "cluster name" : "cluster 0",
                                             "reference mRNA number" : 1,
                                             "alternative mRNA number" : 1}]]
        }

    print("\n*************Testing the 'new' annotation_match function*************")
    
    for test in test_dict:
        
        print(f"\n{test} test")
        
        result = comp.annotation_match(test_dict[test][0], False)
        print(f"result : {result}\n")
        print(f"expected result : {test_dict[test][1]}\n")
        assert result == test_dict[test][1]
        
        
# test function for the 'comparison.py' function 'annotation_match' (main annotation comparison function) with the 'create_strings' parameter as 'True' (uses old program version)
def test_old_annotation_match():
    
    # dictionary of inputs and expected ouputs for each test file for the 'comparison.py' function 'annotation_match' (main annotation comparison function)
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
                 
                 'reverse' : [cl.Cluster(name="cluster 0",
                                         end=299,
                                                    loci={'ref': [locus.Locus(name='chr2A_00611930', 
                                                                             mRNAs={'chr2A_00611930_mrna': [0, 59, 90, 149, 170, 199]}, 
                                                                             start=100, 
                                                                             end=299, 
                                                                             direction='reverse')],
                                                          'alt': [locus.Locus(name='chr2A_00611930', 
                                                                                   mRNAs={'chr2A_00611930_mrna': [0, 59, 90, 149, 170, 199]}, 
                                                                                   start=100, 
                                                                                   end=299, 
                                                                                   direction='reverse')]}),
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
                 
                 'reverse-modified' : [cl.Cluster(name="cluster 0",
                                                  end=299,
                                                    loci={'ref': [locus.Locus(name='chr2A_00611930', 
                                                                             mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                                                             start=100, 
                                                                             end=299, 
                                                                             direction='reverse')],
                                                          'alt': [locus.Locus(name='chr2A_00611930', 
                                                                                   mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 219, 240, 299]}, 
                                                                                   start=100, 
                                                                                   end=299, 
                                                                                   direction='reverse')]}),
                                          [{"reference" : 'chr2A_00611930',
                                             "reference start" : 100,
                                             "reference end" : 299,
                                             "alternative" : 'chr2A_00611930',
                                             "alternative start" : 100,
                                             "alternative end" : 299,
                                             "reference mRNA" : 'chr2A_00611930_mrna',
                                             "alternative mRNA" : 'chr2A_00611930_mrna',
                                             "mismatch/match" : [60, 10, 90],
                                             "identity" : 37.5,
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
        
        result = comp.annotation_match(test_dict[test][0], True, False)
        print(f"result : {result}\n")
        print(f"expected result : {test_dict[test][1]}\n")
        assert result == test_dict[test][1]
        
# return 
# - nb_matches:
# - nb_mismatches_EI:
# - nb_mismatches_RF:
# - list intervals of EI mismatches (upper bounds excluded) 
# - list intervals of RF mismatches (upper bounds included) => inconsistent 
# swith to MatchInfo class      
# test function for the 'comparison.py' function 'compute_matches_mismatches_EI_RF' (CDS intervals comparison function)
def old_test_compute_matches_mismatches_EI_RF():
    
    # dictionary of inputs and expected ouputs for the function
    test_dict = {
        "identical" : [[100, 129, 150, 209, 240, 299],
                       iu.OrderedIntervals(intervals=[100, 129, 150, 209, 240, 299], include_ub=True), 
                       [100, 129, 150, 209, 240, 299],
                       (150, 0, 0, [], [])],
        
        "fusion" : [[100, 129, 150, 209, 240, 299],
                       iu.OrderedIntervals(intervals=[100, 129, 150, 209, 240, 299], include_ub=True), 
                       [100, 209, 240, 299],
                       (30, 20, 120, [130, 150], [150, 209, 240, 299])],
        
        "shift" : [[100, 129, 150, 209, 240, 299],
                       iu.OrderedIntervals(intervals=[100, 129, 150, 209, 240, 299], include_ub=True), 
                       [100, 129, 151, 209, 240, 299],
                       (30, 1, 119, [150, 151], [151, 209, 240, 299])],
        
        "reverse-reverse" : [[0, 59, 90, 149, 170, 199],
                             iu.OrderedIntervals(intervals=[0, 59, 90, 149, 170, 199], include_ub=True), 
                             [0, 59, 90, 149, 170, 199],
                             (150, 0, 0, [], [])],
        
        "diff-start-before" : [[100, 129, 150, 209, 240, 299],
                       iu.OrderedIntervals(intervals=[100, 129, 150, 209, 240, 299], include_ub=True), 
                       [40, 69, 90, 149, 180, 239],
                       (30, 180, 30, [40, 70, 90, 100, 130, 180, 210, 300], [100, 129])],
        
        "diff-start-after" : [[100, 129, 150, 209, 240, 299],
                       iu.OrderedIntervals(intervals=[100, 129, 150, 209, 240, 299], include_ub=True), 
                       [160, 189, 210, 269, 300, 359],
                       (30, 180, 30, [100, 130, 150, 160, 190, 240, 270, 360], [160, 189])],
        
        'overlapping-loci (first cluster/locus)' : [[50, 149],
                       iu.OrderedIntervals(intervals=[50, 149], include_ub=True), 
                       [100, 249],
                       (0, 150, 50, [50, 100, 150, 250], [100, 149])],
        
        'overlapping-loci (first cluster / 2nd locus)' : [[200, 349],
                       iu.OrderedIntervals(intervals=[200, 349], include_ub=True), 
                       [100, 249],
                       (0, 200, 50, [100, 200, 250, 350], [200, 249])],
                         
        "length_computation" : [[8, 13],
                       iu.OrderedIntervals(intervals=[8, 13], include_ub=True), 
                       [1, 4, 12, 13],
                       (2, 8, 0, [1, 5, 8, 12], [])]
        }
    
    print("\n*************Testing the compute_matches_mismatches_EI_RF function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        result = comp.compute_matches_mismatches_EI_RF(test_dict[test][0], test_dict[test][1], test_dict[test][2], False)
        print(f"result : {result}\n")
        print(f"expected result : {test_dict[test][3]}\n")
        assert result == test_dict[test][3]

def test_mismatchInfo_init():
    MI = comp.MismatchInfo([150,150])
    assert MI.nb == 1

def test_compute_matches_mismatches_EI_RF():
    
    # dictionary of inputs and expected outputs for the function
    test_dict = {
        "identical" : [
            [100, 129, 150, 209, 240, 299],
            [100, 129, 150, 209, 240, 299], 
            [150, 0, 0],
            MrnaMatchInfo(matches=150, 
                         mismatches_EI=MismatchInfo([]), 
                         mismatches_RF=MismatchInfo([]), 
                         genomic_overlap=200, 
                         ref_id="ref_mrna_id", 
                         alt_id="alt_mrna_id")
        ],
        
        "fusion" : [
            [100, 129, 150, 209, 240, 299],
            [100, 209, 240, 299],
            [30, 20, 120],
            MrnaMatchInfo(matches=30, 
                         mismatches_EI=MismatchInfo([130, 149]), 
                         mismatches_RF=MismatchInfo([150, 209, 240, 299]), 
                         genomic_overlap=200, 
                         ref_id="ref_mrna_id", 
                         alt_id="alt_mrna_id")
        ],
        
        "shift" : [
            [100, 129, 150, 209, 240, 299],
            [100, 129, 151, 209, 240, 299],
            [30, 1, 119],
            MrnaMatchInfo(matches=30, 
                         mismatches_EI=MismatchInfo([150, 150]), 
                         mismatches_RF=MismatchInfo([151, 209, 240, 299]), 
                         genomic_overlap=200, 
                         ref_id="ref_mrna_id", 
                         alt_id="alt_mrna_id")
        ],
        
        "reverse-reverse" : [
            [0, 59, 90, 149, 170, 199],
            [0, 59, 90, 149, 170, 199],
            [150, 0, 0],
            MrnaMatchInfo(matches=150, 
                         mismatches_EI=MismatchInfo([]), 
                         mismatches_RF=MismatchInfo([]), 
                         genomic_overlap=200, 
                         ref_id="ref_mrna_id", 
                         alt_id="alt_mrna_id")
        ],
        
        "diff-start-before" : [
            [100, 129, 150, 209, 240, 299],
            [40, 69, 90, 149, 180, 239],
            [30, 180, 30],
            MrnaMatchInfo(matches=30, 
                         mismatches_EI=MismatchInfo([40, 69, 90, 99, 130, 179, 210, 299]), 
                         mismatches_RF=MismatchInfo([100, 129]), 
                         genomic_overlap=140, 
                         ref_id="ref_mrna_id", 
                         alt_id="alt_mrna_id")
        ],
        
        "diff-start-after" : [
            [100, 129, 150, 209, 240, 299],
            [160, 189, 210, 269, 300, 359],
            [30, 180, 30],
            MrnaMatchInfo(matches=30, 
                         mismatches_EI=MismatchInfo([100, 129, 150, 159, 190, 239, 270, 359]), 
                         mismatches_RF=MismatchInfo([160, 189]), 
                         genomic_overlap=140, 
                         ref_id="ref_mrna_id", 
                         alt_id="alt_mrna_id")
        ],
        
        "overlapping-loci (first cluster/locus)" : [
            [50, 149], 
            [100, 249],
            [0, 150, 50],
            MrnaMatchInfo(matches=0, 
                         mismatches_EI=MismatchInfo([50, 99, 150, 249]), 
                         mismatches_RF=MismatchInfo([100, 149]), 
                         genomic_overlap=50, 
                         ref_id="ref_mrna_id", 
                         alt_id="alt_mrna_id")
        ],
        
        "overlapping-loci (first cluster / 2nd locus)" : [
            [200, 349], 
            [100, 249],
            [0, 200, 50],
            MrnaMatchInfo(matches=0, 
                         mismatches_EI=MismatchInfo([100, 199, 250, 349]), 
                         mismatches_RF=MismatchInfo([200, 249]), 
                         genomic_overlap=50, 
                         ref_id="ref_mrna_id", 
                         alt_id="alt_mrna_id")
        ],
                         
        "length_computation" : [
            [8, 13], 
            [1, 4, 12, 13],
            [2, 8, 0],
            MrnaMatchInfo(matches=2, 
                         mismatches_EI=MismatchInfo([1, 4, 8, 11]), 
                         mismatches_RF=MismatchInfo([]), 
                         genomic_overlap=6, 
                         ref_id="ref_mrna_id", 
                         alt_id="alt_mrna_id")
        ]
    }
    
    print("\n*************Testing the compute_matches_mismatches_EI_RF function*************")
    
    for test in test_dict:
        print(f"\n{test} file test")
        ref_mrna, alt_mrna, expected_counts, expected_matchInfo = test_dict[test]
        
        # Call the function with proper parameters
        result:MrnaMatchInfo = comp.compute_matches_mismatches_EI_RF("ref_mrna_id", "alt_mrna_id", ref_mrna, alt_mrna, False)
        counts:tuple = [result.matches, result.mismatches_EI.nb, result.mismatches_RF.nb]
        assert counts == expected_counts, f"Counts don't match: {counts} vs {expected_counts}"
        assert result == expected_matchInfo, f"Matches don't match: {result.matches} vs {expected_matchInfo}"



