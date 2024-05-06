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
from main import get_gff_borders, create_vectors, pair_vector_comparison, matrix_to_identity

## This function tests the program on multiple basic 'artificial' test files and checks if their
# return values match what is expected
#
def test_basic():
    
    path1 = "./data/tests/" + "basic_test.gff3"
    path2 = "./data/tests/" + "identical_test.gff3"
    path3 = "./data/tests/" + "minus-CDS_test.gff3"
    path4 = "./data/tests/" + "fusion_test.gff3"
    path5 = "./data/tests/" + "shift_test.gff3"
    path6 = "./data/tests/" + "reverse_test.gff3"
    path7 = "./data/tests/" + "diff-start-before_test.gff3"
    path8 = "./data/tests/" + "diff-start-after_test.gff3"
    
    
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
    
    bord6 = get_gff_borders(path6)
    assert bord6 == {'chr2A_00611930': [300, 240, 210, 150, 130, 100]}
    
    bord7 = get_gff_borders(path7)
    assert bord7 == {'chr2A_00611930': [40, 70, 90, 150, 180, 240]}
    
    bord8 = get_gff_borders(path8)
    assert bord8 == {'chr2A_00611930': [160, 190, 210, 270, 300, 360]}
    
    
    # test of the structure string creation
    
    vect1 = create_vectors( bord1 )
    assert vect1 =={'chr2A_00611930': [100, "12312312312312312312312312312300000000000000000000123123123123123123123123123123123123123123123123123123123123000000000000000000000000000000123123123123123123123123123123123123123123123123123123123123"]}
    
    vect2 = create_vectors( bord2 )
    assert vect2 == {'chr2A_00611930': [100, "12312312312312312312312312312300000000000000000000123123123123123123123123123123123123123123123123123123123123000000000000000000000000000000123123123123123123123123123123123123123123123123123123123123"]}
    
    vect3 = create_vectors( bord3 )
    assert vect3 == {'chr2A_00611930': [100, "12312312312312312312312312312300000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000123123123123123123123123123123123123123123123123123123123123"]}
    
    vect4 = create_vectors( bord4 )
    assert vect4 == {'chr2A_00611930': [100, "12312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312312000000000000000000000000000000312312312312312312312312312312312312312312312312312312312312"]}
    
    vect5 = create_vectors( bord5 )
    assert vect5 == {'chr2A_00611930': [100, "12312312312312312312312312312300000000000000000000012312312312312312312312312312312312312312312312312312312312000000000000000000000000000000312312312312312312312312312312312312312312312312312312312312"]}
    
    vect6 = create_vectors( bord6 )
    assert vect6 =={'chr2A_00611930': [100, "12312312312312312312312312312300000000000000000000123123123123123123123123123123123123123123123123123123123123000000000000000000000000000000123123123123123123123123123123123123123123123123123123123123"]}
    
    vect7 = create_vectors( bord7 )
    assert vect7 =={'chr2A_00611930': [40, "12312312312312312312312312312300000000000000000000123123123123123123123123123123123123123123123123123123123123000000000000000000000000000000123123123123123123123123123123123123123123123123123123123123"]}
    
    vect8 = create_vectors( bord8 )
    assert vect8 =={'chr2A_00611930': [160, "12312312312312312312312312312300000000000000000000123123123123123123123123123123123123123123123123123123123123000000000000000000000000000000123123123123123123123123123123123123123123123123123123123123"]}
    
    
    # test of the vector comparisons of the test files
    
    comp1 = pair_vector_comparison(vect1["chr2A_00611930"], vect2["chr2A_00611930"])
    assert comp1 == [[50, 0, 0, 0],[0, 50, 0, 0],[0, 0, 50, 0],[0, 0, 0, 50]]
    
    comp2 = pair_vector_comparison(vect1["chr2A_00611930"], vect3["chr2A_00611930"])
    assert comp2 == [[50, 0, 0, 0],[20, 30, 0, 0],[20, 0, 30, 0],[20, 0, 0, 30]]
    
    comp3 = pair_vector_comparison(vect1["chr2A_00611930"], vect4["chr2A_00611930"])
    assert comp3 == [[30, 7, 7, 6],[0, 10, 0, 40],[0, 40, 10, 0],[0, 0, 40, 10]]
    
    comp4 = pair_vector_comparison(vect1["chr2A_00611930"], vect5["chr2A_00611930"])
    assert comp4 == [[50, 0, 0, 0],[1, 10, 0, 39],[0, 40, 10, 0],[0, 0, 40, 10]]
    
    comp5 = pair_vector_comparison(vect1["chr2A_00611930"], vect6["chr2A_00611930"])
    assert comp5 == [[50, 0, 0, 0],[0, 50, 0, 0],[0, 0, 50, 0],[0, 0, 0, 50]]
    
    comp6 = pair_vector_comparison(vect1["chr2A_00611930"], vect7["chr2A_00611930"])
    assert comp6 == [[20, 30, 30, 30], [30, 10, 10, 0], [30, 0, 10, 10], [30, 10, 0, 10]]
    
    comp7 = pair_vector_comparison(vect1["chr2A_00611930"], vect8["chr2A_00611930"])
    assert comp7 == [[20, 30, 30, 30], [30, 10, 0, 10], [30, 10, 10, 0], [30, 0, 10, 10]]
    
    
    # test of the identity computation from the pair comparison matrix
    
    ident1 = matrix_to_identity(comp1)
    assert ident1 == 100.0
    
    ident2 = matrix_to_identity(comp2)
    assert ident2 == 60.0
    
    ident3 = matrix_to_identity(comp3)
    assert ident3 == 17.6
    
    ident4 = matrix_to_identity(comp4)
    assert ident4 == 20.0
    
    ident5 = matrix_to_identity(comp5)
    assert ident5 == 100.0
    
    ident6 = matrix_to_identity(comp6)
    assert ident6 == 12.5
    
    ident7 = matrix_to_identity(comp7)
    assert ident7 == 12.5
    
    
