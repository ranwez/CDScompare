#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 11:02:52 2024

@author: vetea
"""

import os, sys

# code adapted from https://csatlas.com/python-import-file-module/#import_a_file_in_a_different_directory
script_dir = os.path.dirname( __file__ )
sys.path.append( script_dir )
import read_files as rf


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
        
        result = rf.get_structure_id(test_dict[test][0], False)
        print(f"result : {result}\n")
        print(test_dict[test][1])
        assert result == test_dict[test][1]


# test function for the 'annot_CSC.py' function 'get_gff_borders' (CDS coordinates acquisition function)
def test_get_gff_borders():
    
    # dictionary of inputs and expected ouputs for each test file for the 'annot_CSC.py' function 'get_gff_borders' (CDS coordinates acquisition function)
    test_dict = {
        "basic" : ["./data/tests/basic_test.gff3", 
                   "chr2A_direct",
                        {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}}],
        
        "identical" : ["./data/tests/identical_test.gff3",
                       "chr2A_direct",
                            {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}}],
        
        "minus-CDS" : ["./data/tests/minus-CDS_test.gff3",
                       "chr2A_direct",
                       {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 129, 240, 299]}}],
        
        "fusion" : ["./data/tests/fusion_test.gff3",
                    "chr2A_direct",
                    {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 209, 240, 299]}}],
        
        "shift" : ["./data/tests/shift_test.gff3",
                   "chr2A_direct",
                   {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 129, 151, 209, 240, 299]}}],
        
        "reverse" : ["./data/tests/reverse_test.gff3",
                     "chr2A_reverse",
                     {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]}}],
        
        "diff-start-before" : ["./data/tests/diff-start-before_test.gff3",
                               "chr2A_direct",
                               {"chr2A_00611930" : {'chr2A_00611930_mrna': [40, 69, 90, 149, 180, 239]}}],
        
        "diff-start-after" : ["./data/tests/diff-start-after_test.gff3",
                              "chr2A_direct",
                              {"chr2A_00611930" : {'chr2A_00611930_mrna': [160, 189, 210, 269, 300, 359]}}],
        
        "multiple mRNAs" : ["./data/tests/multiple_mRNAs_test.gff3",
                            "chr2A_direct",
                            {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 129, 150, 179, 240, 299],
                                                 'chr2A_00611930_mrna.2' : [100, 129, 240, 299]}}],
        
        "basic-2-loci" : ["./data/tests/basic-2-loci_test.gff3",
                          "chr2A_direct",
                              {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]},
                               "chr2A_00620000" : {'chr2A_00620000_mrna': [600, 699, 800, 899]}}],
        
        "identical-2-loci" : ["./data/tests/identical-2-loci_test.gff3",
                              "chr2A_direct",
                              {"chr2A_00611930" : {'chr2A_00611930_mrna': [100, 129, 150, 209, 240, 299]},
                               "chr2A_00620000" : {'chr2A_00620000_mrna': [600, 699, 800, 899]}}],
        
        "overlapping-loci" : ["./data/tests/overlapping-loci_test.gff3",
                              "chr2A_direct",
                              {"chr2A_1000" : {'chr2A_1000_mrna' : [50, 149]},
                               "chr2A_2000" : {'chr2A_2000_mrna' : [200, 349]},
                               "chr2A_3000" : {'chr2A_3000_mrna' : [400, 549]},
                              "chr2A_4000" : {'chr2A_4000_mrna' : [650, 699]},
                              "chr2A_5000" : {'chr2A_5000_mrna' : [750, 799]}}],
        
        "overlapping-loci-alt" : ["./data/tests/overlapping-loci-alt_test.gff3",
                                  "chr2A_direct",
                                  {"chr2A_1000" : {'chr2A_1000_mrna' : [100, 249]},
                                   "chr2A_2000" : {'chr2A_2000_mrna' : [300, 449]},
                                   "chr2A_3000" : {'chr2A_3000_mrna' : [500, 599]},
                                   "chr2A_4000" : {'chr2A_4000_mrna' : [650, 699]},
                                   "chr2A_5000" : {'chr2A_5000_mrna' : [750, 779]},
                                   "chr2A_6000" : {'chr2A_6000_mrna' : [790, 849]}}],
        
        "length_computation_ref" : ["./data/tests/length_computation_ref_test.gff3", 
                                    "chr2A_direct",
                        {"chr2A_00611930" : {'chr2A_00611930_mrna': [8, 13]}}],
        
        "length_computation_alt" : ["./data/tests/length_computation_alt_test.gff3", 
                                    "chr2A_direct",
                        {"chr2A_00611930" : {'chr2A_00611930_mrna': [1, 4, 12, 13]}}]
        
        }
    
    print("\n*************Testing the get_gff_borders function*************")
    
    for test in test_dict:
        print(f"\n{test} file test")
        result = rf.get_gff_borders(test_dict[test][0], False)
        
        # verify presence of expected mRNAs
        dna_mol = test_dict[test][1]
        for loc_id, loc in test_dict[test][2].items():
            print(f"Result keys : {result[dna_mol].keys()}")
            print(f"Expected key : {loc_id}")
            assert loc_id in result[dna_mol]
            expected_mRNA = test_dict[test][2][loc_id]
            print(f"Result mRNAs : {result[dna_mol][loc_id].mRNAs}")
            print(f"Expected mRNA : {test_dict[test][2][loc_id]}")
            assert result[dna_mol][loc_id].contain_mrnas(**expected_mRNA)