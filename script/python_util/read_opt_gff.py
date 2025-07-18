#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 21:04:46 2024

@author: ranwez
"""

import sys
import polars as pl
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

def gff_to_cdsInfo(gff_file: str, relevant_gene_ids: Optional[set[str]] = None) -> dict[str, GeneInfo]:
    """
    Parse a GFF file and return a Polars DataFrame.

    Parameters:
        file_path (str): Path to the GFF file.

    Returns:
        pl.DataFrame: A Polars DataFrame containing the parsed GFF data with added ID, ParentID, mRNA and gene columns.
    """
    column_names = [
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ]

    column_types = {
        "seqid": pl.Categorical,
        "source": pl.Categorical,
        "type": pl.Categorical,
        "start": pl.UInt32, 
        "end": pl.UInt32,
        "score": pl.Utf8,
        "strand": pl.Categorical,
        "phase": pl.UInt8, 
        "attributes": pl.Utf8
    }

    try:
         # 1. get CDS infos
        cds_df = (pl.scan_csv(gff_file, separator="\t", has_header=False, 
                         comment_prefix="#", new_columns=column_names, 
                         schema_overrides=column_types, null_values=".")
              .filter(pl.col("type") == "CDS")
              .with_columns([
                  pl.col("attributes").str.extract(r"Parent=([^;]+)", 1).alias("mrna_id"),
                  pl.col("phase").fill_null(0)

              ])
              .select("mrna_id", "start", "end", "phase")
              .group_by("mrna_id")
                .agg([
                    pl.col("start").alias("start"),
                    pl.col("end").alias("end"),
                    pl.col("phase").alias("phase")
                ])
          .collect())
        
        # 2. get usefull mRNA infos
        mrna_ids_with_cds = cds_df["mrna_id"].unique()
        mrna_df = (pl.scan_csv(gff_file, separator="\t", has_header=False, 
                         comment_prefix="#", new_columns=column_names, 
                         schema_overrides=column_types,null_values=".")
              .filter(pl.col("type") == "mRNA")
              .with_columns([
                  pl.col("attributes").str.extract(r"ID=([^;]+)", 1).alias("mrna_id"),
                  pl.col("attributes").str.extract(r"Parent=([^;]+)", 1).alias("gene_id")
              ]).filter(pl.col("mrna_id").is_in(mrna_ids_with_cds))
              .select("mrna_id", "gene_id", "start", "end")
              .collect())
    
        # 3. get useful gene info for gene we need: seqid start end and phase
        genes_ids_with_cds = mrna_df["gene_id"].unique()
        gene_df = (pl.scan_csv(gff_file, separator="\t", has_header=False, 
                    comment_prefix="#", new_columns=column_names, 
                    schema_overrides=column_types,null_values=".")
              .filter(pl.col("type") == "gene")
              .with_columns([
                  pl.col("attributes").str.extract(r"ID=([^;]+)", 1).alias("gene_id"),
                  (pl.col("seqid").cast(str) + "_" + pl.when(pl.col("strand") == "+").then(pl.lit("direct")).otherwise(pl.lit("reverse"))).alias("chr_strand").cast(pl.Categorical),
                  pl.col("seqid").alias("chr_id")
              ])
              .filter(pl.col("gene_id").is_in(genes_ids_with_cds))
              .select("gene_id", "chr_id", "start", "end", "strand", "chr_strand")
              .sort(["chr_id", "strand", "start", "end"])
              .collect())
    except Exception as e:
        print(f"Error parsing GFF file: {e}")
        exit(1)

        
    # print number of genes
    print(f"Number of genes opt: {gene_df.height}")   
    
    chr_strand_pairs = gene_df.select("chr_strand").unique().to_dicts()
    chrStrand_2_geneInfos = {pair["chr_strand"]: [] for pair in chr_strand_pairs}
    # init geneInfo
    geneId2geneInfos = {}
    mrnaID2geneInfos = {} 
    for row in gene_df.iter_rows(named=True):
        gene_id = row["gene_id"]
        chr_strand =row["chr_strand"]
        if (gene_id in geneId2geneInfos):
            print (f"Error: Gene {gene_id} appears multiple times in the GFF file.")
            sys.exit(1)
        geneInfo =GeneInfo(
            chr=row["chr_id"],
            strand=row["strand"],
            gene_id=gene_id,
            gene_bounds=Bounds(start=row["start"], end=row["end"])
        )
        geneId2geneInfos[gene_id] = geneInfo
        chrStrand_2_geneInfos[chr_strand].append(geneInfo)
    # adds mNRAs to geneInfos 
    mrna_ids = {} 
    for row in mrna_df.iter_rows(named=True):
        # Create MrnaInfo object
        mrna_id=row["mrna_id"]
        gene_id=row["gene_id"]
        if mrna_id in mrna_ids:
            print(f"Error: mRNA {mrna_id} appears multiple times in the GFF file.")
            sys.exit(1)
        if gene_id not in geneId2geneInfos:
            print(f"Error: mRNA {mrna_id} has no associated gene.")
            sys.exit(1)
        mrnaID2geneInfos[mrna_id] = geneId2geneInfos[gene_id]

    for row in cds_df.iter_rows(named=True):
        mrna_id = row["mrna_id"]
        gene = mrnaID2geneInfos.get(mrna_id)     
        cds_bounds = [PhasedBounds(start=start, end=end,phase=phase) for start, end, phase in zip(row["start"], row["end"], row["phase"])]
        cds_bounds.sort(key=lambda x: x.start)  # Sort by start position

        # Vérifier les chevauchements des régions CDS
        for i in range(len(cds_bounds) - 1):
            if cds_bounds[i].end >= cds_bounds[i+1].start:
                print(f"Error: mRNA {mrna_id} has overlapping CDS regions: {cds_bounds[i].end} >= {cds_bounds[i+1].start}")
                exit(1)
        
        # Transformer en liste plate pour utilisation ultérieure
        cds_coords = [coord for bounds in cds_bounds for coord in (bounds.start, bounds.end)]
            
        gene.mRNAs.append(
            MrnaInfo(
                #gene=gene,
                mrna_id=mrna_id,
                cds_bounds=cds_coords,
                phase=row["phase"][0] if gene.strand == "+" else row["phase"][-1],
                mrna_bounds=Bounds(cds_bounds[0].start, cds_bounds[-1].end)
                )
        )
     
    for gene in geneId2geneInfos.values():
        # define the coding bounds using min and max of the mRNA bounds
        if not gene.mRNAs:
            print(f"Error: Gene {gene.gene_id} has no mRNAs.")
            sys.exit(1)
        gene.coding_bounds = Bounds(
            start=min(mrna.mrna_bounds.start for mrna in gene.mRNAs),
            end=max(mrna.mrna_bounds.end for mrna in gene.mRNAs)
        )
    return chrStrand_2_geneInfos