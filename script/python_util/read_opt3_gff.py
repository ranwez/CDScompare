#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 21:04:46 2024

@author: ranwez
"""
STRING_CACHE_DIRECT = "direct" 
STRING_CACHE_REVERSE = "reverse"

import re
from locus import Locus
from array import array
from typing import Optional

def compile_id_regex():
    """Returns compiled regex for extracting ID from GFF attributes"""
    return re.compile(r'\bID=([^;\n\r]+)')

def compile_parent_regex():
    """Returns compiled regex for extracting Parent from GFF attributes"""
    return re.compile(r'\bParent=([^;\n\r]+)')

def gff_to_cdsInfo(gff_file: str, relevant_gene_ids: Optional[set[str]] = None) -> dict[str, Locus]:
    """
    Parse a GFF file and extract gene, mRNA and CDS information, returning Locus objects directly.
    
    Args:
        gff_file: Path to the GFF file
        relevant_gene_ids: Optional set of gene IDs to filter (if None, all genes are processed)
    
    Returns:
        A dictionary mapping chromosome_strand keys to lists of Locus objects
    """
    from locus import Locus
    
    id_regex = compile_id_regex()
    parent_regex = compile_parent_regex()

    mrna_id2CDS = {}      # Dictionary to store CDS by mRNA
    gene_id2mRNA = {}     # Dictionary to store mRNA by gene
    gene_id2info = {}     # Dictionary to store gene information
    
    # Parse the GFF file line by line
    with open(gff_file, "r") as file:
        for line in file:
            if line.startswith("#") or not line.strip():
                continue
            infos = line.strip().split("\t")
            type = infos[2]
            
            
            if type == "CDS":
                start, end = int(infos[3]), int(infos[4])
                mrna_id = parent_regex.search(infos[8]).group(1)
                if mrna_id not in mrna_id2CDS:
                    mrna_id2CDS[mrna_id] = []
                phase = int(infos[7]) if infos[7] != "." else 0
                mrna_id2CDS[mrna_id].append((phase, start, end))
            
            elif type == "mRNA":
                mrna_id = id_regex.search(infos[8]).group(1)
                gene_id = parent_regex.search(infos[8]).group(1)
                
                # Skip if we're only interested in specific genes
                if relevant_gene_ids and gene_id not in relevant_gene_ids:
                    continue
                    
                if gene_id not in gene_id2mRNA:
                    gene_id2mRNA[gene_id] = set()
                gene_id2mRNA[gene_id].add(mrna_id)
            
            elif type == "gene":
                gene_id = id_regex.search(infos[8]).group(1)
                
                # Skip if we're only interested in specific genes
                if relevant_gene_ids and gene_id not in relevant_gene_ids:
                    continue
                # format gene info as a single string
                gene_id2info[gene_id] = (f"{infos[0]}\t{infos[3]}\t{infos[4]}\t{infos[6]}")
   
    print(f"Number of genes (opt3): {len(gene_id2info)}")
    input("Press Enter to continue...")
    
    # Build Locus objects directly
    chrStrand_2_loci = {}
    
    for gene_id, infoStr in gene_id2info.items():
        gene_mrnas = gene_id2mRNA.pop(gene_id, set())  # Récupérer et supprimer immédiatement
        if not gene_mrnas:
            continue
        infos = infoStr.split("\t")
        start, end = int(infos[1]), int(infos[2])
        # Process all mRNAs for this gene
        mrna_dict = {}  # To hold mRNA ID -> CDS bounds
        phases_dict = {}  # To hold mRNA ID -> phase
        
        for mrna_id in gene_mrnas:
            if mrna_id in mrna_id2CDS:
                # Sort CDS by start position
                cds_list = sorted(mrna_id2CDS.pop(mrna_id), key=lambda x: x[1])
                
                # Check for overlapping CDS regions
                for i in range(len(cds_list) - 1):
                    if cds_list[i][2] >= cds_list[i+1][1]:
                        print(f"Error: mRNA {mrna_id} has overlapping CDS regions")
                        exit(1)
                
                # Store CDS bounds and phase
                mrna_dict[mrna_id] = array('L', (coord for cds in cds_list for coord in (cds[1], cds[2])))
                phases_dict[mrna_id] = cds_list[0][0] if infos[3] == "+" else cds_list[-1][0]
        
        # Only create a Locus if we have at least one valid mRNA
        if mrna_dict:
            # Create Locus directly
            locus = Locus(
                name=gene_id,
                start=start,
                end=end,
                mRNAs=mrna_dict,
                phases = phases_dict
            )
            
            # Add to chr_strand dictionary
            direction=STRING_CACHE_DIRECT if infos[3] == "+" else STRING_CACHE_REVERSE
            chr_strand = f"{infos[0]}_{direction}"
            if chr_strand not in chrStrand_2_loci:
                chrStrand_2_loci[chr_strand] = []
            
            chrStrand_2_loci[chr_strand].append(locus)
    
    # Sort loci by position within each chromosome_strand
    for chr_strand in chrStrand_2_loci:
        chrStrand_2_loci[chr_strand].sort(key=lambda l: (l.start, l.end))

    mrna_id2CDS.clear()
    gene_id2mRNA.clear()
    gene_id2info.clear()
    return chrStrand_2_loci