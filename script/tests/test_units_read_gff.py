#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test unitaire pour les fonctions de read_gff.py

@author: ranwez
"""

import sys
import os
import unittest

# Ajouter le chemin du répertoire parent pour pouvoir importer les modules
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
sys.path.insert(0, parent_dir)
python_util_dir = os.path.join(parent_dir, "python_util")
sys.path.insert(0, python_util_dir)

import read_opt2_gff as rg
from read_opt2_gff import GeneInfo, MrnaInfo, Bounds

class TestReadGFF(unittest.TestCase):
    """Test des fonctions de read_gff.py"""
    
    def setUp(self):
        """Initialisation des chemins de fichiers de test"""
        self.data_dir = os.path.join(os.path.dirname(parent_dir), "data", "tests")
        self.basic_test_file = os.path.join(self.data_dir, "basic_test.gff3")
        self.basic_2_loci_test_file = os.path.join(self.data_dir, "basic-2-loci_test.gff3")
    
    def test_parse_gff_basic(self):
        """Test de la fonction parse_gff avec le fichier basic_test.gff3"""
        df = rg.parse_gff(self.basic_test_file)
        
        # Vérifier que le DataFrame a été créé correctement
        self.assertIsNotNone(df)
        
        # Vérifier que les colonnes attendues sont présentes
        expected_columns = ["seqid", "type", "start", "end", 
                           "strand", "phase", "ID", "ParentID", 
                           "mRNA", "gene"]
        for col in expected_columns:
            self.assertIn(col, df.columns)
        
        # Vérifier le nombre de lignes
        self.assertEqual(df.height, 5)  # 1 gene, 1 mRNA, 3 CDS
        
        # Vérifier le contenu de base
        gene_row = df.filter(df["type"] == "gene").to_dicts()[0]
        self.assertEqual(gene_row["ID"], "chr2A_00611930")
        self.assertEqual(gene_row["seqid"], "chr2A")
        
        mrna_row = df.filter(df["type"] == "mRNA").to_dicts()[0]
        self.assertEqual(mrna_row["ID"], "chr2A_00611930_mrna")
        self.assertEqual(mrna_row["ParentID"], "chr2A_00611930")
        
        cds_row = df.filter(df["type"] == "CDS").to_dicts()[0]
        self.assertEqual(cds_row["ParentID"], "chr2A_00611930_mrna")
        self.assertEqual(cds_row["mRNA"], "chr2A_00611930_mrna")
        self.assertEqual(cds_row["gene"], "chr2A_00611930")
    
    def test_gff_to_cdsInfo_basic(self):
        """Test de la fonction gff_to_cdsInfo avec le fichier basic_test.gff3"""
        chrStrand_2_geneInfo = rg.gff_to_cdsInfo(self.basic_test_file)
        
        # Vérifier qu'on a bien un dictionnaire non vide
        self.assertIsInstance(chrStrand_2_geneInfo, dict)
        self.assertEqual(len(chrStrand_2_geneInfo), 1)
        
        # Récupérer l'unique gène
        gene_id = "chr2A_00611930"
        gene= chrStrand_2_geneInfo["chr2A_direct"][0]  
        # Vérifier les propriétés du gène
        self.assertEqual(gene.gene_id, gene_id)
        self.assertEqual(gene.strand, "+")
        self.assertEqual(gene.gene_bounds.start, 100)
        self.assertEqual(gene.gene_bounds.end, 299)
        
        # Vérifier les mRNAs associés
        self.assertEqual(len(gene.mRNAs), 1)
        mrna = gene.mRNAs[0]
        self.assertEqual(mrna.mrna_id, "chr2A_00611930_mrna")
        #self.assertEqual(mrna.gene, gene)  # Vérifier la référence au gène parent
        
        # Vérifier les coordonnées CDS (6 valeurs pour 3 CDS: start1, end1, start2, end2, start3, end3)
        self.assertEqual(len(mrna.cds_bounds), 6)  
        self.assertEqual(mrna.cds_bounds[0], 100)
        self.assertEqual(mrna.cds_bounds[1], 129)
        self.assertEqual(mrna.cds_bounds[2], 150)
        self.assertEqual(mrna.cds_bounds[3], 209)
        self.assertEqual(mrna.cds_bounds[4], 240)
        self.assertEqual(mrna.cds_bounds[5], 299)
    
    def test_gff_to_cdsInfo_2_loci(self):
        """Test de la fonction gff_to_cdsInfo avec le fichier basic-2-loci_test.gff3"""
        chrStrand_2_geneInfo = rg.gff_to_cdsInfo(self.basic_2_loci_test_file)
        
        # Vérifier qu'on a bien un dictionnaire avec 2 gènes
        self.assertIsInstance(chrStrand_2_geneInfo, dict)
        self.assertEqual(len(chrStrand_2_geneInfo), 1)
        self.assertEqual(len(chrStrand_2_geneInfo["chr2A_direct"]), 2)
        
        
        # Vérifier le premier gène
        gene1 = chrStrand_2_geneInfo["chr2A_direct"][0]
        self.assertEqual(gene1.strand, "+")
        self.assertEqual(gene1.gene_bounds.start, 100)
        self.assertEqual(gene1.gene_bounds.end, 299)
        self.assertEqual(len(gene1.mRNAs), 1)
        
        # Vérifier le second gène
        gene2 = chrStrand_2_geneInfo["chr2A_direct"][1] 
        self.assertEqual(gene2.strand, "+")
        self.assertEqual(gene2.gene_bounds.start, 600)
        self.assertEqual(gene2.gene_bounds.end, 899)
        self.assertEqual(len(gene2.mRNAs), 1)
        
        # Vérifier les mRNAs du premier gène
        mrna1 = gene1.mRNAs[0]
        self.assertEqual(mrna1.mrna_id, "chr2A_00611930_mrna")
        
        # Vérifier les mRNAs du second gène
        mrna2 = gene2.mRNAs[0]
        self.assertEqual(mrna2.mrna_id, "chr2A_00620000_mrna")
        
        # Vérifier les coordonnées CDS du premier mRNA
        self.assertEqual(len(mrna1.cds_bounds), 6)  # 3 CDS regions (start, end) x 3
        
        # Vérifier les coordonnées CDS du second mRNA
        self.assertEqual(len(mrna2.cds_bounds), 4)  # 2 CDS regions (start, end) x 2
        self.assertEqual(mrna2.cds_bounds[0], 600)
        self.assertEqual(mrna2.cds_bounds[1], 699)
        self.assertEqual(mrna2.cds_bounds[2], 800)
        self.assertEqual(mrna2.cds_bounds[3], 899)
        
        # Vérifier que les phases sont correctes
        self.assertEqual(mrna1.phase, 0)  # Premier élément pour brin +
        self.assertEqual(mrna2.phase, 0)  # Premier élément pour brin +
    def test_gff_to_cdsInfo_error_multiple_genes(self):
        file="/Users/ranwez/My_data/Projects/2023_ANNOT_WHEAT_LRR/BleTendre/01_sorted_input_gffs/alt_chr_ncbi_sorted.gff"
        dic = rg.gff_to_cdsInfo(file)
        assert ( len(dic) > 1)
        file="/Users/ranwez/My_data/Projects/2023_ANNOT_WHEAT_LRR/BleTendre/01_sorted_input_gffs/ref_chr_urgi_HC_sorted.gff"
        dic = rg.gff_to_cdsInfo(file)
        assert ( len(dic) > 1)
if __name__ == "__main__":
    unittest.main()
