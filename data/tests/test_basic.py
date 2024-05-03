#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:21:54 2024

@author: vetea
"""

import os, sys

# code récupéré/adapté de https://csatlas.com/python-import-file-module/#import_a_file_in_a_different_directory
script_dir = os.path.dirname( __file__ )
main_dir = os.path.join( script_dir, '..', '..', 'script' )
sys.path.append( main_dir )
from main import get_gff_borders, create_vectors, pair_vector_comparison

## This function tests the program on multiple basic 'artificial' test files and checks if their
# return values correspodn to what is expected
#
def test_basic():
    
    path1 = "basic_test.gff3"
    path2 = "identical_test.gff3"
    path3 = "minus-CDS_test.gff3"
    path4 = "fusion_test.gff3"
    path5 = "shift_test.gff3"
    
    # test of the CDS coordinates acquisition
    
    bord1 = get_gff_borders(path1)
    assert bord1 == {'chr2A_00611930': [100, 130, 150, 210, 240, 300]}
    
    bord2 = get_gff_borders(path2)
    assert bord2 == {'chr2A_00611930': [100, 130, 150, 210, 240, 300]}
    
    bord3 = get_gff_borders(path3)
    assert bord3 == {'chr2A_00611930': [100, 130, 240, 300]}
    
    bord4 = get_gff_borders(path4)
    assert bord4 == {'chr2A_00611930': [100, 210, 240, 300]}
    
    bord5 = get_gff_borders(path5)
    assert bord5 == {'chr2A_00611930': [100, 130, 151, 210, 240, 300]}
    
    # test of the structure string creation
    
    vect1 = create_vectors( bord1 )
    assert vect1 =={'chr2A_00611930': "12312312312312312312312312312300000000000000000000123123123123123123123123123123123123123123123123123123123123000000000000000000000000000000123123123123123123123123123123123123123123123123123123123123"}
    
    vect2 = create_vectors( bord2 )
    assert vect2 == {'chr2A_00611930': "12312312312312312312312312312300000000000000000000123123123123123123123123123123123123123123123123123123123123000000000000000000000000000000123123123123123123123123123123123123123123123123123123123123"}
    
    vect3 = create_vectors( bord3 )
    assert vect3 == {'chr2A_00611930': "12312312312312312312312312312300000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000123123123123123123123123123123123123123123123123123123123123"}
    
    vect4 = create_vectors( bord4 )
    assert vect4 == {'chr2A_00611930': "12312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312000000000000000000000000000000312312312312312312312312312312312312312312312312312312312312"}
    
    vect5 = create_vectors( bord5 )
    assert vect5 == {'chr2A_00611930': "12312312312312312312312312312300000000000000000000012312312312312312312312312312312312312312312312312312312312000000000000000000000000000000312312312312312312312312312312312312312312312312312312312312"}
    
    # test of the vector comparisons of the test files
    
    assert(pair_vector_comparison(vect1["chr2A_00611930"], vect2["chr2A_00611930"])) == [[50, 0, 0, 0],[0, 50, 0, 0],[0, 0, 50, 0],[0, 0, 0, 50]]
    
    assert(pair_vector_comparison(vect1["chr2A_00611930"], vect3["chr2A_00611930"])) == [[50, 0, 0, 0],[20, 30, 0, 0],[20, 0, 30, 0],[20, 0, 0, 30]]
    
    assert(pair_vector_comparison(vect1["chr2A_00611930"], vect4["chr2A_00611930"])) == [[30, 7, 7, 6],[0, 10, 0, 40],[0, 40, 10, 0],[0, 0, 40, 10]]
    
    assert(pair_vector_comparison(vect1["chr2A_00611930"], vect5["chr2A_00611930"])) == [[50, 0, 0, 0],[1, 10, 0, 39],[0, 40, 10, 0],[0, 0, 40, 10]]