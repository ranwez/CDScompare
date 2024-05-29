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

## This script tests the 'annot_CSC.py' program on multiple basic 'artificial' test files and checks if their return values match what is expected


# test function for the 'annot_CSC.py' function 'get_structure_id' (structure id acquisition function)
def test_get_structure_id():
    test_dict = {        
        "TRITD_HC" : ["chr2A	PGSB_Jan2017	gene	22128	24635	.	-	.	ID=TRITD2Av1G000030;primconf=HC\n",
                      "TRITD2Av1G000030"],
        
        "EXP_chr2A" : ["chr2A	exonerate:protein2genome:local	gene	608664	611579	.	-	.	ID=chr2A_00611930;color=2;comment=modified:yes / Gene-Class:Non-canonical / start_changed_to3f;note=Origin:OSJnip_Chr04_04489784 / pred:prot2genome / prot-%25-ident:51.6 / prot-%25-cov:99.2118 / exo_corr:NA / Origin-Fam:NLR / Origin-Class:Canonical / noStart / Gene-Class:Non-canonical\n",
                       "chr2A_00611930"],
        
        "chr2A" : ["chr2A	exonerate:protein2genome:local	gene	608664	611630	.	-	.	ID=chr2A_00611930;comment=Origin:OSJnip_Chr04_04489784 / pred:prot2genome / prot-%-ident:51.6 / prot-%-cov:99.2118 / exo_corr:NA / Origin-Fam:NLR / Origin-Class:Canonical / noStart / Gene-Class:Non-canonical;color=2\n",
                   "chr2A_00611930"],
        
        "annot_best" : ["contig	exonerate:protein2genome:local	gene	5131	9328	.	+	.	ID=DWSvevo3_contig_0000005131;comment=Origin:DWSvevo1_chr2A_4125497 / pred:prot2genome / prot-%-ident:44.063 / prot-%-cov:68.2268 / score:51.7466 / scoreNC:46.9466 / exo_corr:corrected_intron / exo_corr:modif_stop / exo_corr:fix_overlap / Origin-Fam:Non-canonical / Origin-Class:LRR-RLK Gene-Class:Non-canonical / noStart / pbFrameshift / unexpectedSplicingSite / stopInFrame;color=2\n",
                        "DWSvevo3_contig_0000005131"],
        
        "basic_test" : ["chr2A	exonerate:protein2genome:local	gene	100	300	.	-	.	ID=chr2A_00611930\n",
                        "chr2A_00611930"]
        }
    
    print("\n*************Testing the get_structure_id function*************")
    
    for test in test_dict:
        
        print(f"\n{test} test")
        
        result = annot_CSC.get_structure_id(test_dict[test][0], False, False)
        print(result)
        print(test_dict[test][1])
        assert result == test_dict[test][1]


# test function for the 'annot_CSC.py' function 'get_gff_borders' (CDS coordinates acquisition function)
def test_get_gff_borders():
    
    # dictionary of inputs and expected ouputs for each test file for the 'annot_CSC.py' function 'get_gff_borders' (CDS coordinates acquisition function)
    test_dict = {
        "basic" : ["./data/tests/basic_test.gff3", 
                        {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]}}],
        
        "identical" : ["./data/tests/identical_test.gff3",
                            {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]}}],
        
        "minus-CDS" : ["./data/tests/minus-CDS_test.gff3",
                       {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 130, 240, 300]}}],
        
        "fusion" : ["./data/tests/fusion_test.gff3",
                    {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 210, 240, 300]}}],
        
        "shift" : ["./data/tests/shift_test.gff3",
                   {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 130, 151, 210, 240, 300]}}],
        
        "reverse" : ["./data/tests/reverse_test.gff3",
                     {"chr2A_00611930" : {'chr2A_00611930_mrna': [300, 240, 210, 150, 130, 100]}}],
        
        "diff-start-before" : ["./data/tests/diff-start-before_test.gff3",
                               {"chr2A_00611930" : {'chr2A_00611930_mrna': [40, 70, 90, 150, 180, 240]}}],
        
        "diff-start-after" : ["./data/tests/diff-start-after_test.gff3",
                              {"chr2A_00611930" : {'chr2A_00611930_mrna': [160, 190, 210, 270, 300, 360]}}],
        
        "basic-2-loci" : ["./data/tests/basic-2-loci_test.gff3",
                              {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]},
                               "chr2A_00620000" : {'chr2A_00620000_mrna': [600, 700, 800, 900]}}],
        
        "identical-2-loci" : ["./data/tests/identical-2-loci_test.gff3",
                              {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]},
                               "chr2A_00620000" : {'chr2A_00620000_mrna': [600, 700, 800, 900]}}],
        
        "overlapping-loci" : ["./data/tests/overlapping-loci_test.gff3",
                              {"chr2A_1000" : {'chr2A_1000_mrna' : [50,150]},
                               "chr2A_2000" : {'chr2A_2000_mrna' : [200, 350]},
                               "chr2A_3000" : {'chr2A_3000_mrna' : [400, 550]},
                              "chr2A_4000" : {'chr2A_4000_mrna' : [650, 700]},
                              "chr2A_5000" : {'chr2A_5000_mrna' : [750, 800]}}],
        
        "overlapping-loci-alt" : ["./data/tests/overlapping-loci-alt_test.gff3",
                                  {"chr2A_1000" : {'chr2A_1000_mrna' : [100,250]},
                                   "chr2A_2000" : {'chr2A_2000_mrna' : [300, 450]},
                                   "chr2A_3000" : {'chr2A_3000_mrna' : [500, 600]},
                                   "chr2A_4000" : {'chr2A_4000_mrna' : [650, 700]},
                                   "chr2A_5000" : {'chr2A_5000_mrna' : [750, 780]},
                                   "chr2A_6000" : {'chr2A_6000_mrna' : [790, 850]}}],
        }
    
    print("\n*************Testing the get_gff_borders function*************")
    
    for test in test_dict:
        print(f"\n{test} file test")
        result = annot_CSC.get_gff_borders(test_dict[test][0], False, False)
        
        for loc_id, loc in test_dict[test][1].items():
            print(f"Result keys : {result.keys()}")
            print(f"Expected key : {loc_id}")
            assert loc_id in result
            expected_mRNA = test_dict[test][1][loc_id]
            print(f"Result mRNAs : {result[loc_id].mRNAs}")
            print(f"Expected mRNA : {test_dict[test][1][loc_id]}")
            assert result[loc_id].contain_mrnas(**expected_mRNA)
        
        
# test function for the 'annot_CSC.py' function 'annotation_sort' (locus list creation and sorting function)
def test_annotation_sort():
    
    # dictionary of inputs and expected ouputs for each test file for the 'annot_CSC.py' function 'annotation_sort' (locus list creation and sorting function)
    test_dict = {"overlapping-loci" : [
                                {'chr2A_1000': annot_CSC.Locus(name='chr2A_1000', 
                                                     mRNAs={'chr2A_1000_mrna': [50, 150]}, 
                                                     start=50, 
                                                     end=150, 
                                                     direction='direct'), 
                                'chr2A_2000': annot_CSC.Locus(name='chr2A_2000', 
                                                     mRNAs={'chr2A_2000_mrna': [200, 350]}, 
                                                     start=200, 
                                                     end=350, 
                                                     direction='direct'), 
                                'chr2A_3000': annot_CSC.Locus(name='chr2A_3000', 
                                                     mRNAs={'chr2A_3000_mrna': [400, 550]}, 
                                                     start=400, 
                                                     end=550, 
                                                     direction='direct'),
                                'chr2A_4000': annot_CSC.Locus(name='chr2A_4000', 
                                                     mRNAs={'chr2A_4000_mrna': [650, 700]}, 
                                                     start=650, 
                                                     end=700, 
                                                     direction='direct'),
                                'chr2A_5000': annot_CSC.Locus(name='chr2A_5000', 
                                                     mRNAs={'chr2A_5000_mrna': [750, 800]}, 
                                                     start=750, 
                                                     end=800, 
                                                     direction='direct')},
                                                    
                                {'chr2A_1000': annot_CSC.Locus(name='chr2A_1000', 
                                                     mRNAs={'chr2A_1000_mrna': [100, 250]}, 
                                                     start=100, 
                                                     end=250, 
                                                     direction='direct'), 
                                'chr2A_2000': annot_CSC.Locus(name='chr2A_2000', 
                                                     mRNAs={'chr2A_2000_mrna': [300, 450]}, 
                                                     start=300, 
                                                     end=450, 
                                                     direction='direct'), 
                                'chr2A_3000': annot_CSC.Locus(name='chr2A_3000', 
                                                     mRNAs={'chr2A_3000_mrna': [500, 600]}, 
                                                     start=500, 
                                                     end=600, 
                                                     direction='direct'),
                                'chr2A_4000': annot_CSC.Locus(name='chr2A_4000', 
                                                     mRNAs={'chr2A_4000_mrna': [650, 700]}, 
                                                     start=650, 
                                                     end=700, 
                                                     direction='direct'),
                                'chr2A_5000': annot_CSC.Locus(name='chr2A_5000', 
                                                     mRNAs={'chr2A_5000_mrna': [750, 780]}, 
                                                     start=750, 
                                                     end=780, 
                                                     direction='direct'),
                                'chr2A_6000': annot_CSC.Locus(name='chr2A_6000', 
                                                     mRNAs={'chr2A_6000_mrna': [790, 850]}, 
                                                     start=790, 
                                                     end=850, 
                                                     direction='direct')},
                                                      
                                [(50, 150, 'chr2A_1000', True), 
                               (100, 250, 'chr2A_1000', False), 
                               (200, 350, 'chr2A_2000', True), 
                               (300, 450, 'chr2A_2000', False), 
                               (400, 550, 'chr2A_3000', True), 
                               (500, 600, 'chr2A_3000', False), 
                               (650, 700, 'chr2A_4000', False), 
                               (650, 700, 'chr2A_4000', True), 
                               (750, 780, 'chr2A_5000', False), 
                               (750, 800, 'chr2A_5000', True), 
                               (790, 850, 'chr2A_6000', False)]]
    }
    
    print("\n*************Testing the annotation_sort function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        result = annot_CSC.annotation_sort(test_dict[test][0], test_dict[test][1], False, False)
        print(result)
        print(test_dict[test][2])
        assert result == test_dict[test][2]
        
        
# test function for the 'annot_CSC.py' function 'construct_clusters' (overlapping loci grouping function)
def test_construct_clusters():
    
    # dictionary of inputs and expected ouputs for each test file for the 'annot_CSC.py' function 'construct_clusters' (overlapping loci grouping function)
    test_dict = {'simple' : [{
        'chr2A_1000': annot_CSC.Locus(name='chr2A_1000', 
                                 mRNAs={'chr2A_1000_mrna': [50, 150]}, 
                                 start=50, 
                                 end=150, 
                                 direction='direct'), 
        'chr2A_2000': annot_CSC.Locus(name='chr2A_2000', 
                                 mRNAs={'chr2A_2000_mrna': [200, 350]}, 
                                 start=200, 
                                 end=350, 
                                 direction='direct'), 
        'chr2A_3000': annot_CSC.Locus(name='chr2A_3000', 
                                 mRNAs={'chr2A_3000_mrna': [400, 550]}, 
                                 start=400, 
                                 end=550, 
                                 direction='direct'),
        'chr2A_4000': annot_CSC.Locus(name='chr2A_4000', 
                                 mRNAs={'chr2A_4000_mrna': [650, 700]}, 
                                 start=650, 
                                 end=700, 
                                 direction='direct'),
        'chr2A_5000': annot_CSC.Locus(name='chr2A_5000', 
                                 mRNAs={'chr2A_5000_mrna': [750, 800]}, 
                                 start=750, 
                                 end=800, 
                                 direction='direct')},
                    
        {'chr2A_1000': annot_CSC.Locus(name='chr2A_1000', 
                             mRNAs={'chr2A_1000_mrna': [100, 250]}, 
                             start=100, 
                             end=250, 
                             direction='direct'), 
        'chr2A_2000': annot_CSC.Locus(name='chr2A_2000', 
                             mRNAs={'chr2A_2000_mrna': [300, 450]}, 
                             start=300, 
                             end=450, 
                             direction='direct'), 
        'chr2A_3000': annot_CSC.Locus(name='chr2A_3000', 
                             mRNAs={'chr2A_3000_mrna': [500, 600]}, 
                             start=500, 
                             end=600, 
                             direction='direct'),
        'chr2A_4000': annot_CSC.Locus(name='chr2A_4000', 
                             mRNAs={'chr2A_4000_mrna': [650, 700]}, 
                             start=650, 
                             end=700, 
                             direction='direct'),
        'chr2A_5000': annot_CSC.Locus(name='chr2A_5000', 
                             mRNAs={'chr2A_5000_mrna': [750, 780]}, 
                             start=750, 
                             end=780, 
                             direction='direct'),
        'chr2A_6000': annot_CSC.Locus(name='chr2A_6000', 
                             mRNAs={'chr2A_6000_mrna': [790, 850]}, 
                             start=790, 
                             end=850, 
                             direction='direct')},
        
        [(50, 150, 'chr2A_1000', True), 
          (100, 250, 'chr2A_1000', False), 
          (200, 350, 'chr2A_2000', True), 
          (300, 450, 'chr2A_2000', False), 
          (400, 550, 'chr2A_3000', True), 
          (500, 600, 'chr2A_3000', False), 
          (650, 700, 'chr2A_4000', False), 
          (650, 700, 'chr2A_4000', True), 
          (750, 780, 'chr2A_5000', False), 
          (750, 800, 'chr2A_5000', True), 
          (790, 850, 'chr2A_6000', False)],
        
        {'cluster 0':  [{'chr2A_1000_mrna': [50, 150]}, {'chr2A_2000_mrna': [200, 350]}, {'chr2A_3000_mrna': [400, 550]}], 
          'cluster 6':  [{'chr2A_4000_mrna': [650, 700]}], 
          'cluster 8':  [{'chr2A_5000_mrna': [750, 800]}]},
        
        {'cluster 0':  [{'chr2A_1000_mrna': [100, 250]}, {'chr2A_2000_mrna': [300, 450]}, {'chr2A_3000_mrna': [500, 600]}], 
          'cluster 6':  [{'chr2A_4000_mrna': [650, 700]}], 
          'cluster 8':  [{'chr2A_5000_mrna': [750, 780]}, {'chr2A_6000_mrna': [790, 850]}]}]
        }

    print("\n*************Testing the construct_clusters function*************")
    
    for test in test_dict:
        
        print(f"\n{test} test")
        
        result = annot_CSC.construct_clusters(test_dict[test][0], test_dict[test][1], test_dict[test][2], False, False)
        for cluster_id, cluster in test_dict[test][3].items():
            for mRNA in cluster:
                result_mRNAs = result.get_mRNAs(cluster_id, "ref")
                print(f"Expected mRNA = {mRNA}")
                print(f"Result mRNAs = {result_mRNAs}")
                assert mRNA in result_mRNAs
        for cluster_id, cluster in test_dict[test][4].items():
            for mRNA in cluster:
                result_mRNAs = result.get_mRNAs(cluster_id, "alt")
                print(f"Expected mRNA = {mRNA}")
                print(f"Result mRNAs = {result_mRNAs}")
                assert mRNA in result_mRNAs
    

# test function for the 'annot_CSC.py' function 'get_area_bounds' (CDS coordinates fusion function)
def test_get_area_bounds():
    
    # dictionary of inputs and expected ouputs for each test file for the 'annot_CSC.py' function 'get_area_bounds' (CDS coordinates fusion function)
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
        
        result = annot_CSC.get_area_bounds(test_dict[test][0], test_dict[test][1], False, False)
        print(result)
        print(test_dict[test][2])
        assert result == test_dict[test][2]
    

# test function for the 'annot_CSC.py' function 'is_in_cds' (CDS inclusion in comparison area function)
def test_is_in_cds():
    
    # dictionary of inputs and expected ouputs for each test file for the 'annot_CSC.py' function 'is_in_cds' (CDS inclusion in comparison area function)
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
        
        result = annot_CSC.is_in_cds(test_dict[test][0], test_dict[test][1], False, False)
        print(result)
        print(test_dict[test][2])
        assert result == test_dict[test][2]


# test function for the 'annot_CSC.py' function 'compare_loci' (new locus comparison function)
def test_compare_loci():
    
    # dictionary of inputs and expected ouputs for each test file for the 'annot_CSC.py' function 'compare_loci' (new locus comparison function)
    test_dict = {
        "identical" : [annot_CSC.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]}, 
                                 start=100, 
                                 end=300, 
                                 direction='direct'), 
                       annot_CSC.Locus(name='chr2A_00611930', 
                                    mRNAs={'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]}, 
                                    start=100, 
                                    end=300, 
                                    direction='direct'),
                       ([0, 150], 100.0, [])],
        
        "minus-CDS" : [annot_CSC.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]}, 
                                 start=100, 
                                 end=300, 
                                 direction='direct'),
                       annot_CSC.Locus(name='chr2A_00611930', 
                                    mRNAs={'chr2A_00611930_mrna': [100, 130, 240, 300]}, 
                                    start=100, 
                                    end=300, 
                                    direction='direct'),
                       ([60, 90], 60.0, ['150-210'])],
        
        "fusion" : [annot_CSC.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]}, 
                                 start=100, 
                                 end=300, 
                                 direction='direct'),
                    annot_CSC.Locus(name='chr2A_00611930', 
                                mRNAs={'chr2A_00611930_mrna': [100, 210, 240, 300]}, 
                                start=100, 
                                end=300, 
                                direction='direct'),
                    ([140, 30], 17.6, ['130-150', '150-210', '240-300'])],
        
        "shift" : [annot_CSC.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]}, 
                                 start=100, 
                                 end=300, 
                                 direction='direct'),
                   annot_CSC.Locus(name='chr2A_00611930', 
                                mRNAs={'chr2A_00611930_mrna': [100, 130, 151, 210, 240, 300]}, 
                                start=100, 
                                end=300, 
                                direction='direct'),
                   ([120, 30], 20.0, ['150-151', '151-210', '240-300'])],
        
        "reverse-reverse" : [annot_CSC.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [300, 240, 210, 150, 130, 100]}, 
                                 start=300, 
                                 end=100, 
                                 direction='reverse'),
                     annot_CSC.Locus(name='chr2A_00611930', 
                                mRNAs={'chr2A_00611930_mrna': [300, 240, 210, 150, 130, 100]}, 
                                start=300, 
                                end=100, 
                                direction='reverse'),
                     ([0, 150], 100.0, [])],
        
        "diff-start-before" : [annot_CSC.Locus(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]}, 
                                            start=100, 
                                            end=300, 
                                            direction='direct'),
                               annot_CSC.Locus(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [40, 70, 90, 150, 180, 240]}, 
                                            start=40, 
                                            end=240, 
                                            direction='direct'),
                               ([210, 30], 12.5, ['40-70', '90-100', '100-130', '130-150', '150-180', '210-240', '240-300'])],
        
        "diff-start-after" : [annot_CSC.Locus(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]}, 
                                            start=100, 
                                            end=300, 
                                            direction='direct'),
                              annot_CSC.Locus(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [160, 190, 210, 270, 300, 360]}, 
                                            start=160, 
                                            end=360, 
                                            direction='direct'),
                              ([210, 30], 12.5, ['100-130', '150-160', '160-190', '190-210', '210-240', '270-300', '300-360'])],
        
        "basic-2-loci (second locus)" : [annot_CSC.Locus(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [600, 700, 800, 900]}, 
                                            start=600, 
                                            end=900, 
                                            direction='direct'),
                                         annot_CSC.Locus(name='chr2A_00611930', 
                                                    mRNAs={'chr2A_00611930_mrna': [600, 700, 800, 900]}, 
                                                    start=600, 
                                                    end=900, 
                                                    direction='direct'),
                                         ([0, 200], 100.0, [])]
        }
    
    print("\n*************Testing the compare_loci function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        result = annot_CSC.compare_loci(test_dict[test][0], test_dict[test][1], False, False)
        print(f"Expected result = {test_dict[test][2]}")
        print(f"Result = {result}")
        assert result == test_dict[test][2]

    
# test function for the 'annot_CSC.py' function 'create_vectors' (structure string creation function)
def test_create_vectors():
    
    # dictionary of inputs and expected ouputs for each test file for the 'annot_CSC.py' function 'create_vectors' (structure string creation function)
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
        
        result = annot_CSC.create_vectors(test_dict[test][0], False, False)
        print(result)
        print(test_dict[test][1])
        assert result == test_dict[test][1]
     
        
# test function for the 'annot_CSC.py' function 'create_vectors' (structure string creation function)
def test_old_compare_loci():
    
    # dictionary of inputs and expected ouputs for each test file for the 'annot_CSC.py' function 'old_compare_loci' (old locus comparison function)
    test_dict = {
        "identical" : [annot_CSC.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]}, 
                                 start=100, 
                                 end=300, 
                                 direction='direct'), 
                       annot_CSC.Locus(name='chr2A_00611930', 
                                    mRNAs={'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]}, 
                                    start=100, 
                                    end=300, 
                                    direction='direct'),
                       ([0, 150], 100.0)],
        
        "minus-CDS" : [annot_CSC.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]}, 
                                 start=100, 
                                 end=300, 
                                 direction='direct'),
                       annot_CSC.Locus(name='chr2A_00611930', 
                                    mRNAs={'chr2A_00611930_mrna': [100, 130, 240, 300]}, 
                                    start=100, 
                                    end=300, 
                                    direction='direct'),
                       ([60, 90], 60.0)],
        
        "fusion" : [annot_CSC.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]}, 
                                 start=100, 
                                 end=300, 
                                 direction='direct'),
                    annot_CSC.Locus(name='chr2A_00611930', 
                                mRNAs={'chr2A_00611930_mrna': [100, 210, 240, 300]}, 
                                start=100, 
                                end=300, 
                                direction='direct'),
                    ([140, 30], 17.6)],
        
        "shift" : [annot_CSC.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]}, 
                                 start=100, 
                                 end=300, 
                                 direction='direct'),
                   annot_CSC.Locus(name='chr2A_00611930', 
                                mRNAs={'chr2A_00611930_mrna': [100, 130, 151, 210, 240, 300]}, 
                                start=100, 
                                end=300, 
                                direction='direct'),
                   ([120, 30], 20.0)],
        
        "reverse-reverse" : [annot_CSC.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [300, 240, 210, 150, 130, 100]}, 
                                 start=300, 
                                 end=100, 
                                 direction='reverse'),
                     annot_CSC.Locus(name='chr2A_00611930', 
                                mRNAs={'chr2A_00611930_mrna': [300, 240, 210, 150, 130, 100]}, 
                                start=300, 
                                end=100, 
                                direction='reverse'),
                     ([0, 150], 100.0)],
        
        "diff-start-before" : [annot_CSC.Locus(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]}, 
                                            start=100, 
                                            end=300, 
                                            direction='direct'),
                               annot_CSC.Locus(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [40, 70, 90, 150, 180, 240]}, 
                                            start=40, 
                                            end=240, 
                                            direction='direct'),
                               ([210, 30], 12.5)],
        
        "diff-start-after" : [annot_CSC.Locus(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]}, 
                                            start=100, 
                                            end=300, 
                                            direction='direct'),
                              annot_CSC.Locus(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [160, 190, 210, 270, 300, 360]}, 
                                            start=160, 
                                            end=360, 
                                            direction='direct'),
                              ([210, 30], 12.5)],
        
        "basic-2-loci (second locus)" : [annot_CSC.Locus(name='chr2A_00611930', 
                                            mRNAs={'chr2A_00611930_mrna': [600, 700, 800, 900]}, 
                                            start=600, 
                                            end=900, 
                                            direction='direct'),
                                         annot_CSC.Locus(name='chr2A_00611930', 
                                                    mRNAs={'chr2A_00611930_mrna': [600, 700, 800, 900]}, 
                                                    start=600, 
                                                    end=900, 
                                                    direction='direct'),
                                         ([0, 200], 100.0)],
        
        'overlapping-loci (first cluster/locus)' : [annot_CSC.Locus(name='chr2A_1000', 
                                                             mRNAs={'chr2A_1000_mrna': [50, 150]}, 
                                                             start=50, 
                                                             end=150, 
                                                             direction='direct'),
                                                    annot_CSC.Locus(name='chr2A_1000', 
                                                              mRNAs={'chr2A_1000_mrna': [100, 250]}, 
                                                              start=100, 
                                                              end=250, 
                                                              direction='direct'),
                                                    ([200, 0], 0.0)],
        
        'overlapping-loci (first cluster / 2nd locus)' : [annot_CSC.Locus(name='chr2A_2000', 
                                                             mRNAs={'chr2A_2000_mrna': [200, 350]}, 
                                                             start=200, 
                                                             end=350, 
                                                             direction='direct'),
                                                          annot_CSC.Locus(name='chr2A_2000', 
                                                              mRNAs={'chr2A_2000_mrna': [300, 450]}, 
                                                              start=300, 
                                                              end=450, 
                                                              direction='direct'),
                                                          ([250, 0], 0.0)],
        }
    
    print("\n*************Testing the old_compare_loci function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        result = annot_CSC.old_compare_loci(test_dict[test][0], test_dict[test][1], False, False)
        print(result)
        print(test_dict[test][2])
        assert result == test_dict[test][2]

        
# test function for the 'annot_CSC.py' function 'annotation_match' (main annotation comparison function) with the 'create_strings' parameter as 'False' (uses new program version)
def test_new_annotation_match():
    
    # dictionary of inputs and expected ouputs for each test file for the 'annot_CSC.py' function 'annotation_match' (main annotation comparison function)
    test_dict = {'basic' : [[annot_CSC.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]}, 
                                 start=100, 
                                 end=300, 
                                 direction='direct')],
                             [annot_CSC.Locus(name='chr2A_00611930', 
                                                mRNAs={'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]}, 
                                                start=100, 
                                                end=300, 
                                                direction='direct')],
                             [{"reference" : 'chr2A_00611930',
                                "reference start" : 100,
                                "reference end" : 300,
                                "alternative" : 'chr2A_00611930',
                                "alternative start" : 100,
                                "alternative end" : 300,
                                "mismatch/match" : [0, 150],
                                "identity" : 100.0,
                                "mismatch zones" : []}]],
                 
                 'shift' : [[annot_CSC.Locus(name='chr2A_00611930', 
                                              mRNAs={'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]}, 
                                              start=100, 
                                              end=300, 
                                              direction='direct')],
                            [annot_CSC.Locus(name='chr2A_00611930', 
                                                mRNAs={'chr2A_00611930_mrna': [100, 130, 151, 210, 240, 300]}, 
                                                start=100, 
                                                end=300, 
                                                direction='direct')],
                            [{"reference" : 'chr2A_00611930',
                            "reference start" : 100,
                            "reference end" : 300,
                            "alternative" : 'chr2A_00611930',
                            "alternative start" : 100,
                            "alternative end" : 300,
                            "mismatch/match" : [120, 30],
                            "identity" : 20.0,
                            "mismatch zones" : ['150-151', '151-210', '240-300']}]],
                 
                 'overlapping-loci (first cluster)' : [[annot_CSC.Locus(name='chr2A_1000', 
                                              mRNAs={'chr2A_1000_mrna': [50, 150]}, 
                                              start=50, 
                                              end=150, 
                                              direction='direct'), 
                                        annot_CSC.Locus(name='chr2A_2000', 
                                              mRNAs={'chr2A_2000_mrna': [200, 350]}, 
                                              start=200, 
                                              end=350, 
                                              direction='direct'), 
                                        annot_CSC.Locus(name='chr2A_3000', 
                                              mRNAs={'chr2A_3000_mrna': [400, 550]}, 
                                              start=400, 
                                              end=550, 
                                              direction='direct')],
                                 
                                       [annot_CSC.Locus(name='chr2A_1000', 
                                          mRNAs={'chr2A_1000_mrna': [100, 250]}, 
                                          start=100, 
                                          end=250, 
                                          direction='direct'), 
                                        annot_CSC.Locus(name='chr2A_2000', 
                                          mRNAs={'chr2A_2000_mrna': [300, 450]}, 
                                          start=300, 
                                          end=450, 
                                          direction='direct'), 
                                        annot_CSC.Locus(name='chr2A_3000', 
                                          mRNAs={'chr2A_3000_mrna': [500, 600]}, 
                                          start=500, 
                                          end=600, 
                                          direction='direct')],
                                    [{"reference" : 'chr2A_2000',
                                    "reference start" : 200,
                                    "reference end" : 350,
                                    "alternative" : 'chr2A_2000',
                                    "alternative start" : 300,
                                    "alternative end" : 450,
                                    "mismatch/match" : [250, 0],
                                    "identity" : 0.0,
                                    "mismatch zones" : ['200-300', '300-350', '350-450']},
                                     {"reference" : 'chr2A_1000',
                                    "reference start" : 50,
                                    "reference end" : 150,
                                    "alternative" : 'chr2A_1000',
                                    "alternative start" : 100,
                                    "alternative end" : 250,
                                    "mismatch/match" : [200, 0],
                                    "identity" : 0.0,
                                    "mismatch zones" : ['50-100', '100-150', '150-250']},
                                     {"reference" : 'chr2A_3000',
                                     "reference start" : 400,
                                     "reference end" : 550,
                                     "alternative" : 'chr2A_3000',
                                     "alternative start" : 500,
                                     "alternative end" : 600,
                                     "mismatch/match" : [200, 0],
                                     "identity" : 0.0,
                                     "mismatch zones" : ['400-500', '500-550', '550-600']}]]
        }

    print("\n*************Testing the annotation_match function*************")
    
    for test in test_dict:
        
        print(f"\n{test} test")
        
        result = annot_CSC.annotation_match(test_dict[test][0], test_dict[test][1], False, False, False)
        print(result)
        print(test_dict[test][2])
        assert result == test_dict[test][2]
        
        
# test function for the 'annot_CSC.py' function 'annotation_match' (main annotation comparison function) with the 'create_strings' parameter as 'True' (uses old program version)
def test_old_annotation_match():
    
    # dictionary of inputs and expected ouputs for each test file for the 'annot_CSC.py' function 'annotation_match' (main annotation comparison function)
    test_dict = {'basic' : [[annot_CSC.Locus(name='chr2A_00611930', 
                                 mRNAs={'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]}, 
                                 start=100, 
                                 end=300, 
                                 direction='direct')],
                             [annot_CSC.Locus(name='chr2A_00611930', 
                                                mRNAs={'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]}, 
                                                start=100, 
                                                end=300, 
                                                direction='direct')],
                             [{"reference" : 'chr2A_00611930',
                                "reference start" : 100,
                                "reference end" : 300,
                                "alternative" : 'chr2A_00611930',
                                "alternative start" : 100,
                                "alternative end" : 300,
                                "mismatch/match" : [0, 150],
                                "identity" : 100.0,
                                "mismatch zones" : '?'}]],
                 
                 'shift' : [[annot_CSC.Locus(name='chr2A_00611930', 
                                              mRNAs={'chr2A_00611930_mrna': [100, 130, 150, 210, 240, 300]}, 
                                              start=100, 
                                              end=300, 
                                              direction='direct')],
                            [annot_CSC.Locus(name='chr2A_00611930', 
                                                mRNAs={'chr2A_00611930_mrna': [100, 130, 151, 210, 240, 300]}, 
                                                start=100, 
                                                end=300, 
                                                direction='direct')],
                            [{"reference" : 'chr2A_00611930',
                            "reference start" : 100,
                            "reference end" : 300,
                            "alternative" : 'chr2A_00611930',
                            "alternative start" : 100,
                            "alternative end" : 300,
                            "mismatch/match" : [120, 30],
                            "identity" : 20.0,
                            "mismatch zones" : '?'}]],
                 
                 'overlapping-loci (first cluster)' : [[annot_CSC.Locus(name='chr2A_1000', 
                                              mRNAs={'chr2A_1000_mrna': [50, 150]}, 
                                              start=50, 
                                              end=150, 
                                              direction='direct'), 
                                        annot_CSC.Locus(name='chr2A_2000', 
                                              mRNAs={'chr2A_2000_mrna': [200, 350]}, 
                                              start=200, 
                                              end=350, 
                                              direction='direct'), 
                                        annot_CSC.Locus(name='chr2A_3000', 
                                              mRNAs={'chr2A_3000_mrna': [400, 550]}, 
                                              start=400, 
                                              end=550, 
                                              direction='direct')],
                                 
                                       [annot_CSC.Locus(name='chr2A_1000', 
                                          mRNAs={'chr2A_1000_mrna': [100, 250]}, 
                                          start=100, 
                                          end=250, 
                                          direction='direct'), 
                                        annot_CSC.Locus(name='chr2A_2000', 
                                          mRNAs={'chr2A_2000_mrna': [300, 450]}, 
                                          start=300, 
                                          end=450, 
                                          direction='direct'), 
                                        annot_CSC.Locus(name='chr2A_3000', 
                                          mRNAs={'chr2A_3000_mrna': [500, 600]}, 
                                          start=500, 
                                          end=600, 
                                          direction='direct')],
                                    [{"reference" : 'chr2A_2000',
                                    "reference start" : 200,
                                    "reference end" : 350,
                                    "alternative" : 'chr2A_2000',
                                    "alternative start" : 300,
                                    "alternative end" : 450,
                                    "mismatch/match" : [250, 0],
                                    "identity" : 0.0,
                                    "mismatch zones" : '?'},
                                     {"reference" : 'chr2A_1000',
                                    "reference start" : 50,
                                    "reference end" : 150,
                                    "alternative" : 'chr2A_1000',
                                    "alternative start" : 100,
                                    "alternative end" : 250,
                                    "mismatch/match" : [200, 0],
                                    "identity" : 0.0,
                                    "mismatch zones" : '?'},
                                     {"reference" : 'chr2A_3000',
                                     "reference start" : 400,
                                     "reference end" : 550,
                                     "alternative" : 'chr2A_3000',
                                     "alternative start" : 500,
                                     "alternative end" : 600,
                                     "mismatch/match" : [200, 0],
                                     "identity" : 0.0,
                                     "mismatch zones" : '?'}]]
        }

    print("\n*************Testing the annotation_match function*************")
    
    for test in test_dict:
        
        print(f"\n{test} test")
        
        result = annot_CSC.annotation_match(test_dict[test][0], test_dict[test][1], True, False, False)
        print(result)
        print(test_dict[test][2])
        assert result == test_dict[test][2]






