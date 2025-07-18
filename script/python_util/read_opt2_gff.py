#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 21:04:46 2024

@author: ranwez
"""

import sys
import re
from typing import Optional
from attrs import define, field


@define(slots=True)
class Bounds:
    start: int
    end: int

    def __init__(self, start, end):
        self.start=min(start,end)
        self.end=max(start,end)

    def length(self) -> int:
        return self.end - self.start + 1

    def overlap(self, other: "Bounds") -> int:
        if self.end < other.start or other.end < self.start:
            return 0
        return min(self.end, other.end) - max(self.start, other.start) + 1

@define(slots=True, frozen=True)
class PhasedBounds:
    start: int
    end: int
    phase: int = 0


@define(slots=True)
class MrnaInfo:
    #gene: "GeneInfo"
    mrna_id: str
    cds_bounds: list[int] = field(factory=list)
    phase: int = 0
    mrna_bounds: Bounds = None
     

@define(slots=True)
class GeneInfo:
    strand: str
    chr: str
    gene_id: str
    gene_bounds: Bounds
    coding_bounds: Bounds = None
    mRNAs: list[MrnaInfo] = field(factory=list)

    def into_locus(self):
        """
        Convert GeneInfo to Locus object.
        
        Returns:
            Locus: A Locus object containing the gene information.
        """
        from locus import Locus
        return Locus(
            name=self.gene_id,
            mRNAs={mrna.mrna_id: mrna.cds_bounds for mrna in self.mRNAs},
            start=self.gene_bounds.start,
            end=self.gene_bounds.end,
            direction= "direct" if self.strand == "+" else "reverse",
            #phases=[mrna.phase for mrna in self.mRNAs]
            phases={mrna.mrna_id: mrna.phase for mrna in self.mRNAs},
        )

def compile_id_regex():
    """Returns compiled regex for extracting ID from GFF attributes"""
    return re.compile(r'\bID=([^;\n\r]+)')

def compile_parent_regex():
    """Returns compiled regex for extracting Parent from GFF attributes"""
    return re.compile(r'\bParent=([^;\n\r]+)')

def gff_to_cdsInfo(gff_file: str, relevant_gene_ids: Optional[set[str]] = None) -> dict[str, GeneInfo]:
    """
    Parse a GFF file and extract gene, mRNA and CDS information.
    
    Args:
        gff_file: Path to the GFF file
        relevant_gene_ids: Optional set of gene IDs to filter (if None, all genes are processed)
    
    Returns:
        A dictionary mapping chromosome_strand keys to lists of GeneInfo objects
    """
    id_regex = compile_id_regex()
    parent_regex = compile_parent_regex()

    mrna_id2CDS = {}      # Dictionary to store CDS by mRNA
    gene_id2mRNA = {}     # Dictionary to store mRNA by gene
    gene_id2geneinfos = {}  # Dictionary to store GeneInfo objects
    
    # Parse the GFF file line by line
    with open(gff_file, "r") as file:
        for line in file:
            if line.startswith("#") or not line.strip():
                continue
            infos = line.strip().split("\t")
            type = infos[2]
            
            if type in ("CDS", "mRNA", "gene"):
                start = min(int(infos[3]), int(infos[4]))
                end = max(int(infos[3]), int(infos[4]))
            
            if type == "CDS":
                mrna_id = parent_regex.search(infos[8]).group(1)
                if mrna_id not in mrna_id2CDS:
                    mrna_id2CDS[mrna_id] = []
                phase = int(infos[7]) if infos[7] != "." else 0
                mrna_id2CDS[mrna_id].append((phase, start, end))
            
            elif type == "mRNA":
                mrna_id = id_regex.search(infos[8]).group(1)
                gene_id = parent_regex.search(infos[8]).group(1)
                if gene_id not in gene_id2mRNA:
                    gene_id2mRNA[gene_id] = []
                gene_id2mRNA[gene_id].append((mrna_id, start, end))
            
            elif type == "gene":
                gene_id = id_regex.search(infos[8]).group(1)
                gene_id2geneinfos[gene_id] = { "chr": infos[0], "strand": infos[6], "start": start, "end": end }
    
    # Filter genes without mRNAs
    genes_with_mrna = set(gene_id2mRNA.keys())
    genes_to_remove = set(gene_id2geneinfos.keys()) - genes_with_mrna
    for gene_id in genes_to_remove:
        gene_id2geneinfos.pop(gene_id)
    
    # Filter mRNAs without CDS
    mrnas_with_cds = set(mrna_id2CDS.keys())
    
    # For each gene, only keep mRNAs with CDS
    for gene_id, mrnas in list(gene_id2mRNA.items()):
        # Filter mRNAs without CDS for this gene
        gene_id2mRNA[gene_id] = [
            mrna for mrna in mrnas if mrna[0] in mrnas_with_cds
        ]
        
        # If after filtering there are no mRNAs left, remove the gene
        if not gene_id2mRNA[gene_id]:
            if gene_id in gene_id2geneinfos:
                gene_id2geneinfos.pop(gene_id)
    # print number of genes with mRNAs
    print(f"Number of genes with mRNAs (opt2): {len(gene_id2geneinfos)}")

    # Build GeneInfo objects from collected data
    chrStrand_2_geneInfos = {}
    
    for gene_id, gene_data in gene_id2geneinfos.items():
        # Create GeneInfo object
        gene_info = GeneInfo(
            chr=gene_data["chr"],
            strand=gene_data["strand"],
            gene_id=gene_id,
            gene_bounds=Bounds(start=gene_data["start"], end=gene_data["end"])
        )
        
        # Add mRNAs to gene
        for mrna_tuple in gene_id2mRNA.get(gene_id, []):
            mrna_id = mrna_tuple[0]
            if mrna_id in mrna_id2CDS:
                # Sort CDS by start position
                cds_list = sorted(mrna_id2CDS[mrna_id], key=lambda x: x[1])
                
                # Check for overlapping CDS regions
                for i in range(len(cds_list) - 1):
                    if cds_list[i][2] >= cds_list[i+1][1]:
                        print(f"Error: mRNA {mrna_id} has overlapping CDS regions")
                        exit(1)
                
                # Create flat list of CDS coordinates
                cds_bounds = []
                for phase, cds_start, cds_end in cds_list:
                    cds_bounds.extend([cds_start, cds_end])
                
                # Add mRNA to the gene's mRNA list
                mrna_info = MrnaInfo(
                    mrna_id=mrna_id,
                    cds_bounds=cds_bounds,
                    phase=cds_list[0][0] if gene_data["strand"] == "+" else cds_list[-1][0],
                    mrna_bounds=Bounds(cds_list[0][1], cds_list[-1][2])
                )
                gene_info.mRNAs.append(mrna_info)
        
        # Calculate coding bounds for the gene
        if gene_info.mRNAs:
            gene_info.coding_bounds = Bounds(
                start=min(mrna.mrna_bounds.start for mrna in gene_info.mRNAs),
                end=max(mrna.mrna_bounds.end for mrna in gene_info.mRNAs)
            )
        
        # Add gene to the dictionary indexed by chromosome and strand
        chr_strand = f"{gene_data['chr']}_{'direct' if gene_data['strand'] == '+' else 'reverse'}"
        if chr_strand not in chrStrand_2_geneInfos:
            chrStrand_2_geneInfos[chr_strand] = []
        
        chrStrand_2_geneInfos[chr_strand].append(gene_info)
    
    # Sort genes by position within each chromosome_strand
    for chr_strand in chrStrand_2_geneInfos:
        chrStrand_2_geneInfos[chr_strand].sort(key=lambda g: (g.gene_bounds.start, g.gene_bounds.end))
    
    return chrStrand_2_geneInfos