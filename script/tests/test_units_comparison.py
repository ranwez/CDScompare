# -*- coding: utf-8 -*-

import os, sys
import script.python_util.comparison as comp
import script.python_util.cluster as cl
from script.python_util.comparison import MrnaMatchInfo, MismatchInfo
from script.python_util.locus import Locus

# helper function to convert a list of CDS bounds into a Locus object
def cds2locus(cds_list, name):
    """Convert a list of CDS into a Locus object."""
    mrna_dict = {}
    mrna_dict[name] = cds_list
    start = min(cds_list[0],cds_list[-1])
    end = max(cds_list[0],cds_list[-1])
    return Locus.builder(name=name, mRNAs=mrna_dict, start=start, end=end)

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
        "identical" : [Locus.builder(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                 start=100, 
                                 end=299, 
                                 direction='direct'), 
                       Locus.builder(name='chr2A_00611930', 
                                    mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                    start=100, 
                                    end=299, 
                                    direction='direct'),
                        MrnaMatchInfo(matches=150, 
                            mismatches_EI=MismatchInfo([]), 
                            mismatches_RF=MismatchInfo([]), 
                            genomic_overlap=200, 
                            ref_id="chr2A_00611930_mrna", 
                            alt_id="chr2A_00611930_mrna")],
        
        "minus-CDS" : [Locus.builder(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                 start=100, 
                                 end=299, 
                                 direction='direct'),
                       Locus.builder(name='chr2A_00611930', 
                                    mRNAs={'chr2A_00611930_mrna': [100, 129, 240, 299]}, 
                                    start=100, 
                                    end=299, 
                                    direction='direct'),
                          MrnaMatchInfo(matches=90,
                            mismatches_EI=MismatchInfo([150, 209]), 
                            mismatches_RF=MismatchInfo([]), 
                            genomic_overlap=200,
                            ref_id="chr2A_00611930_mrna",
                            alt_id="chr2A_00611930_mrna")],
        
        "fusion" : [Locus.builder(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                 start=100, 
                                 end=299, 
                                 direction='direct'),
                    Locus.builder(name='chr2A_00611930', 
                                mRNAs={'chr2A_00611930_mrna': [100, 209, 240, 299]}, 
                                start=100, 
                                end=299, 
                                direction='direct'),
                    MrnaMatchInfo(matches=30,
                            mismatches_EI=MismatchInfo([130, 149]),
                            mismatches_RF=MismatchInfo([150, 209, 240, 299]),
                            genomic_overlap=200,
                            ref_id="chr2A_00611930_mrna",
                            alt_id="chr2A_00611930_mrna")],
        
        "shift" : [Locus.builder(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                 start=100, 
                                 end=299, 
                                 direction='direct'),
                   Locus.builder(name='chr2A_00611930', 
                                mRNAs={'chr2A_00611930_mrna': [100, 129, 151, 209, 240, 299]}, 
                                start=100, 
                                end=299, 
                                direction='direct'),
                     MrnaMatchInfo(matches=30,
                            mismatches_EI=MismatchInfo([150, 150]),
                            mismatches_RF=MismatchInfo([151, 209, 240, 299]),
                            genomic_overlap=200,
                            ref_id="chr2A_00611930_mrna",
                            alt_id="chr2A_00611930_mrna")],
        
        "reverse" : [Locus.builder(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [0, 59, 90, 149, 170, 199]}, 
                                 start=100, 
                                 end=299, 
                                 direction='reverse'),
                     Locus.builder(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [0, 59, 90, 149, 170, 199]}, 
                                 start=100, 
                                 end=299, 
                                 direction='reverse'),
                     MrnaMatchInfo(matches=150,
                            mismatches_EI=MismatchInfo([]),
                            mismatches_RF=MismatchInfo([]),
                            genomic_overlap=200,
                            ref_id="chr2A_00611930_mrna",
                            alt_id="chr2A_00611930_mrna")],                      
        
        "reverse-modified" : [Locus.builder(name='chr2A_00611930', 
                                          mRNAs={'chr2A_00611930_mrna': [0, 59, 90, 149, 170, 199]}, 
                                          start=100, 
                                          end=299, 
                                          direction='reverse'),
                              Locus.builder(name='chr2A_00611930', 
                                          mRNAs={'chr2A_00611930_mrna': [0, 59, 80, 149, 170, 199]}, 
                                          start=100, 
                                          end=299, 
                                          direction='reverse'),
                              MrnaMatchInfo(matches=90,
                            mismatches_EI=MismatchInfo([80, 89]),
                            mismatches_RF=MismatchInfo([0, 59]),
                            genomic_overlap=200,
                            ref_id="chr2A_00611930_mrna",
                            alt_id="chr2A_00611930_mrna")],
        
        "reverse-basic" : [Locus.builder(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                 start=100, 
                                 end=299, 
                                 direction='direct'),
                            Locus.builder(name='chr2A_00611930', 
                                mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                start=100, 
                                end=299, 
                                direction='reverse'),
                     MrnaMatchInfo(matches=150,
                            mismatches_EI=MismatchInfo([]),
                            mismatches_RF=MismatchInfo([]),
                            genomic_overlap=200,
                            ref_id="chr2A_00611930_mrna",
                            alt_id="chr2A_00611930_mrna")],

        "diff-start-before" : [Locus.builder(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                            start=100, 
                                            end=299, 
                                            direction='direct'),
                               Locus.builder(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [40, 69, 90, 149, 180, 239]}, 
                                            start=40, 
                                            end=239, 
                                            direction='direct'),
                                 MrnaMatchInfo(matches=30,
                            mismatches_EI=MismatchInfo([40, 69, 90, 99, 130, 179, 210, 299]),
                            mismatches_RF=MismatchInfo([100, 129]),
                            genomic_overlap=140,
                            ref_id="chr2A_00611930_mrna",
                            alt_id="chr2A_00611930_mrna")],
        
        "diff-start-after" : [Locus.builder(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                            start=100, 
                                            end=299, 
                                            direction='direct'),
                              Locus.builder(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [160, 189, 210, 269, 300, 359]}, 
                                            start=160, 
                                            end=359, 
                                            direction='direct'),
                                MrnaMatchInfo(matches=30,
                                    mismatches_EI=MismatchInfo([100, 129, 150, 159, 190, 239, 270, 359]),
                                    mismatches_RF=MismatchInfo([160, 189]),
                                    genomic_overlap=140,
                                    ref_id="chr2A_00611930_mrna",
                                    alt_id="chr2A_00611930_mrna")],
        
        "multiple_mRNAs" : [Locus.builder(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                 start=100, 
                                 end=299, 
                                 direction='direct'), 
                       Locus.builder(name='chr2A_00611930', 
                                    mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 179, 240, 299],
                                           'chr2A_00611930_mrna.2': [100, 129, 240, 299]}, 
                                    start=100, 
                                    end=299, 
                                    direction='direct'),
                          MrnaMatchInfo(matches=120,
                                mismatches_EI=MismatchInfo([180, 209]),
                                mismatches_RF=MismatchInfo([]),
                                genomic_overlap=200,
                                ref_id="chr2A_00611930_mrna",
                                alt_id="chr2A_00611930_mrna")],
        
        "basic-2-loci (second locus)" : [Locus.builder(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [600, 699, 800, 899]}, 
                                            start=600, 
                                            end=899, 
                                            direction='direct'),
                                         Locus.builder(name='chr2A_00611930', 
                                                    mRNAs={
                                                        'chr2A_00611930_mrna': [600, 699],
                                                        'chr2A_00611930_mrna.2': [600, 699, 800, 899]
                                                        }, 
                                                    start=600, 
                                                    end=899, 
                                                    direction='direct'),
                                            MrnaMatchInfo(matches=200,
                                                mismatches_EI=MismatchInfo([]),
                                                mismatches_RF=MismatchInfo([]),
                                                genomic_overlap=300,
                                                ref_id="chr2A_00611930_mrna",
                                                alt_id="chr2A_00611930_mrna.2")],
        
        "length_computation" : [Locus.builder(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [8, 13]}, 
                                 start=8, 
                                 end=13, 
                                 direction='direct'), 
                       Locus.builder(name='chr2A_00611930', 
                                    mRNAs={'chr2A_00611930_mrna': [1, 4, 12, 13]}, 
                                    start=1, 
                                    end=13, 
                                    direction='direct'),
                            MrnaMatchInfo(matches=2,
                                mismatches_EI=MismatchInfo([1, 4, 8, 11]),
                                mismatches_RF=MismatchInfo([]),
                                genomic_overlap=6,
                                ref_id="chr2A_00611930_mrna",
                                alt_id="chr2A_00611930_mrna")],
        }
    
    print("\n*************Testing the compare_loci function*************")
    
    for test in test_dict:
        print(f"\n{test} test")
        ref_locus = test_dict[test][0]
        alt_locus = test_dict[test][1]
        cluster_end = max(ref_locus.end, alt_locus.end)+1
        cluster = cl.Cluster(name="cluster 0", loci_ref=[ref_locus], loci_alt=[alt_locus], end=cluster_end)
        reversed = test.find('reverse') != -1
        if (reversed):
            cluster.reverse_loci_coord()
        cluster.set_intervals()
        result = comp.compare_loci(ref_locus, alt_locus)
        if (reversed):
            (rev_mismatches_EI, rev_mismatches_RF) = comp.reverse_coord((result.mismatches_EI.zones, result.mismatches_RF.zones), cluster_end)
            result = MrnaMatchInfo(matches=result.matches,
                                mismatches_EI=MismatchInfo(zones=rev_mismatches_EI),
                                mismatches_RF=MismatchInfo(zones=rev_mismatches_RF),
                                genomic_overlap=result.genomic_overlap,
                                ref_id=result.ref_id,
                                alt_id=result.alt_id)
        print(f"Expected result = {test_dict[test][2]}")
        print(f"Result = {result}")
        assert result == test_dict[test][2]
        
           
        
# test function for the 'comparison.py' function 'annotation_match' (main annotation comparison function) with the 'create_strings' parameter as 'False' (uses new program version)
def test_new_annotation_match():
    
    # dictionary of inputs and expected ouputs for each test file for the 'comparison.py' function 'annotation_match' (main annotation comparison function)
    test_dict = {'basic' : [cl.Cluster(name="cluster 0",
                                       loci_ref= [Locus.builder(name='chr2A_00611930', 
                                                                    mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                                                    start=100, 
                                                                    end=299, 
                                                                    direction='direct')],
                                        loci_alt= [Locus.builder(name='chr2A_00611930', 
                                                                          mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                                                          start=100, 
                                                                          end=299, 
                                                                          direction='direct')]),
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
                                                    loci_ref= [Locus.builder(name='chr2A_00611930', 
                                                                                 mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                                                                 start=100, 
                                                                                 end=299, 
                                                                                 direction='direct')],
                                                    loci_alt= [Locus.builder(name='chr2A_00611930', 
                                                                                       mRNAs={'chr2A_00611930_mrna': [100, 129, 151, 209, 240, 299]}, 
                                                                                       start=100, 
                                                                                       end=299, 
                                                                                       direction='direct')]),
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
                            "mismatch zones" : ([150, 150], [151, 209, 240, 299]),
                            "cluster name" : "cluster 0",
                            "reference mRNA number" : 1,
                            "alternative mRNA number" : 1}]],
                 
                 'reverse' : [cl.Cluster(name="cluster 0",
                                         end=299,
                                                    loci_ref= [Locus.builder(name='chr2A_00611930', 
                                                                             mRNAs={'chr2A_00611930_mrna': [0, 59, 90, 149, 170, 199]}, 
                                                                             start=100, 
                                                                             end=299, 
                                                                             direction='reverse')],
                                                          loci_alt= [Locus.builder(name='chr2A_00611930', 
                                                                                   mRNAs={'chr2A_00611930_mrna': [0, 59, 90, 149, 170, 199]}, 
                                                                                   start=100, 
                                                                                   end=299, 
                                                                                   direction='reverse')]),
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
                                            loci_ref= [Locus.builder(name='chr2A_00611930', 
                                                                     mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                                                     start=100, 
                                                                     end=299, 
                                                                     direction='reverse')],
                                                  loci_alt= [Locus.builder(name='chr2A_00611930', 
                                                                           mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 219, 240, 299]}, 
                                                                           start=100, 
                                                                           end=299, 
                                                                           direction='reverse')]),
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
                                             "mismatch zones" : ([210, 219], [100, 129, 150, 209]),
                                             "cluster name" : "cluster 0",
                                             "reference mRNA number" : 1,
                                             "alternative mRNA number" : 1}]],
                 
                 'multiple_mRNAs' : [cl.Cluster(name="cluster 0",
                                                loci_ref= [Locus.builder(name='chr2A_00611930', 
                                                                         mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}, 
                                                                         start=100, 
                                                                         end=299, 
                                                                         direction='direct')],
                                                      loci_alt= [Locus.builder(name='chr2A_00611930', 
                                                                         mRNAs={'chr2A_00611930_mrna': [100, 129, 150, 179, 240, 299],
                                                                                'chr2A_00611930_mrna.2': [100, 129, 240, 299]}, 
                                                                         start=100, 
                                                                         end=299, 
                                                                         direction='direct')]),
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
                                             "mismatch zones" : ([180, 209], []),
                                             "cluster name" : "cluster 0",
                                             "reference mRNA number" : 1,
                                             "alternative mRNA number" : 2}]],
                 
                 'overlapping-loci (first cluster)' : [cl.Cluster(name="cluster 0",
                                                                  loci_ref=[Locus.builder(name='chr2A_1000', 
                                              mRNAs={'chr2A_1000_mrna': [50, 149]}, 
                                              start=50, 
                                              end=149, 
                                              direction='direct'), 
                                        Locus.builder(name='chr2A_2000', 
                                              mRNAs={'chr2A_2000_mrna': [200, 349]}, 
                                              start=200, 
                                              end=349, 
                                              direction='direct'), 
                                        Locus.builder(name='chr2A_3000', 
                                              mRNAs={'chr2A_3000_mrna': [400, 549]}, 
                                              start=400, 
                                              end=549, 
                                              direction='direct')],
                                                                        loci_alt=[Locus.builder(name='chr2A_1000', 
                                          mRNAs={'chr2A_1000_mrna': [100, 249]}, 
                                          start=100, 
                                          end=249, 
                                          direction='direct'), 
                                        Locus.builder(name='chr2A_2000', 
                                          mRNAs={'chr2A_2000_mrna': [300, 449]}, 
                                          start=300, 
                                          end=449, 
                                          direction='direct'), 
                                        Locus.builder(name='chr2A_3000', 
                                          mRNAs={'chr2A_3000_mrna': [500, 599]}, 
                                          start=500, 
                                          end=599, 
                                          direction='direct')]),
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
                                   "mismatch zones" : ([50, 99, 150, 249], [100, 149]),
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
                                    "mismatch zones" : ([200, 299, 350, 449], [300, 349]),
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
                                     "mismatch zones" : ([400, 499, 550, 599], [500, 549]),
                                     "cluster name" : "cluster 0",
                                     "reference mRNA number" : 1,
                                     "alternative mRNA number" : 1}]],
                 
                 'length_computation' : [cl.Cluster(name="cluster 0",
                                                    loci_ref=[Locus.builder(name='chr2A_00611930', 
                                          mRNAs={'chr2A_00611930_mrna': [8, 13]}, 
                                          start=8, 
                                          end=13, 
                                          direction='direct')], 
                                                          loci_alt=[Locus.builder(name='chr2A_00611930', 
                                             mRNAs={'chr2A_00611930_mrna': [1, 4, 12, 13]}, 
                                             start=1, 
                                             end=13, 
                                             direction='direct')]),
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
                                             "mismatch zones" : ([1, 4, 8, 11], []),
                                             "cluster name" : "cluster 0",
                                             "reference mRNA number" : 1,
                                             "alternative mRNA number" : 1}]]
        }

    print("\n*************Testing the 'new' annotation_match function*************")
    
    for test in test_dict:
        
        print(f"\n{test} test")
        is_reversed= test.find('reverse') != -1
        result = comp.annotation_match(test_dict[test][0], is_reversed,True)
        print(f"result : {result}\n")
        print(f"expected result : {test_dict[test][1]}\n")
        assert result == test_dict[test][1]
        
        
        

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
        ref_locus = cds2locus( ref_mrna, "ref_mrna_id")
        alt_locus = cds2locus( alt_mrna, "alt_mrna_id")
        reversed = test.find('reverse') != -1
        cluster_end = max(ref_locus.end, alt_locus.end)+1
        if (reversed):
            ref_locus.reverse(cluster_end)
            alt_locus.reverse(cluster_end)
        ref_locus.set_mrna_intervals()
        alt_locus.set_mrna_intervals()
        result:MrnaMatchInfo = comp.compute_matches_mismatches_EI_RF(0, 0, ref_locus, alt_locus)
        counts:tuple = [result.matches, result.mismatches_EI.nb, result.mismatches_RF.nb]
        assert counts == expected_counts, f"Counts don't match: {counts} vs {expected_counts}"
        #assert result == expected_matchInfo, f"Matches don't match: {result} vs {expected_matchInfo}"



