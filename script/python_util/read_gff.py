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
    gene: "GeneInfo"
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

def parse_gff(file_path):
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
        "seqid": pl.Utf8,
        "source": pl.Utf8,
        "type": pl.Utf8,
        "start": pl.Int64,
        "end": pl.Int64,
        "score": pl.Utf8,
        "strand": pl.Utf8,
        "phase": pl.Int64,
        "attributes": pl.Utf8,
    }

    try:
        # Lecture lazy avec scan_csv
        lazy_df = pl.scan_csv(
            file_path,
            separator="\t",
            has_header=False,
            comment_prefix="#",
            new_columns=column_names,
            schema_overrides=column_types,
            null_values=".",
        )
        
        # Filtrer en mode lazy
        lazy_df = lazy_df.filter(pl.col("type").is_in(["gene", "mRNA", "CDS"]))

                # Extract IDs from attributes column - AVANT collect()
        # Correction: utiliser pl.col() au lieu de lazy_df["attributes"]
        lazy_df = lazy_df.with_columns([
            pl.col("attributes").str.extract(r"ID=([^;]+)", 1).alias("ID"),
            pl.col("attributes").str.extract(r"Parent=([^;]+)", 1).alias("ParentID"),
        ])

        # Populate mRNA column - ÉGALEMENT AVANT collect()
        lazy_df = lazy_df.with_columns([
            pl.when(pl.col("type") == "mRNA")
            .then(pl.col("ID"))
            .otherwise(pl.col("ParentID"))
            .alias("mRNA")
        ])

        # Create gene column based on mRNA column
        # For gene rows, gene is ID; for others, we need to map mRNA to gene
        lazy_df = lazy_df.with_columns([
            pl.when(pl.col("type") == "gene")
            .then(pl.col("ID"))
            .otherwise(None)
            .alias("gene_temp")
        ])
        
        # Garder uniquement les colonnes nécessaires pour la suite du traitement
        lazy_df = lazy_df.select(
            "type", "start", "end", "strand", "phase", "seqid",
        "ID", "ParentID", "mRNA", "gene_temp"
        )
        df = lazy_df.collect()
        
        # Create mapping for mRNA to gene relationships
        mRNA_mapping = (
            df.filter(df["type"] == "mRNA")
            .select(["ID", "ParentID"])
            .rename({"ID": "mRNA", "ParentID": "gene"})
        )
        
        # Join to populate gene column based on mRNA
        df = df.join(mRNA_mapping, on="mRNA", how="left")

        # Assign gene IDs to gene features
        df = df.with_columns(
            [
                pl.when(df["type"] == "gene")
                .then(df["ID"])
                .otherwise(df["gene"])
                .alias("gene")
            ]
        )
        return df
    except Exception as e:
        print(f"Error parsing GFF file: {e}")
        exit(1)

def gff_to_cdsInfo(gff_file: str, relevant_gene_ids: Optional[set[str]] = None) -> dict[str, GeneInfo]:
    """
    Parse a GFF file and return a dictionary of GeneInfo objects.
    
    Parameters:
        gff_file (str): Path to the GFF file
        relevant_gene_ids (Optional[set[str]]): Set of gene IDs to filter on, or None
    
    Returns:
        dict[str, GeneInfo]: Dictionary of GeneInfo objects keyed by gene ID
    """
    df = parse_gff(gff_file)
    print(f"Parsing GFF file: {gff_file}")
    
    # Filter to keep only feature of relevant genes
    if relevant_gene_ids is not None:
        df = df.filter(df["gene"].is_in(list(relevant_gene_ids)))

    # we focus on first CDS frame, in this case if no phase is given 0 is the correct default
    df = df.with_columns(pl.col("phase").fill_null(0))
        
    # Get CDS coordinates
    cds_coordinates = (
        df.filter(df["type"] == "CDS")
        #.sort("start")
        .group_by("mRNA")
        .agg([
            pl.col("start").alias("start"), 
            pl.col("end").alias("end"),
            pl.col("phase").alias("phase"),
            pl.col("seqid").first().alias("chr_id"),
            pl.col("strand").first().alias("strand"),
            pl.col("gene").first().alias("gene")
        ])
    )
    
    gene_coordinates = (
        df.filter(df["type"] == "gene")
        .group_by("gene")
        .agg([
            pl.col("start").first().alias("gene_start"),
            pl.col("end").first().alias("gene_end"),
            pl.col("strand").first().alias("strand"),
            pl.col("seqid").first().alias("chr_id")
        ])
    )
    
    del df
    # Join CDS and gene information
    merged_data = cds_coordinates.join(gene_coordinates, on="gene", how="left")
    del cds_coordinates
    del gene_coordinates
    print(f"Number of genes: {merged_data.select(pl.col('gene')).unique().height}")
    
    merged_data = merged_data.with_columns(
        (pl.col("chr_id") + "_" + pl.when(pl.col("strand") == "+").then(pl.lit("direct")).otherwise(pl.lit("reverse"))).alias("chr_strand")
    )
    # get the pairs of chromosome and strand
    chr_strand_pairs = merged_data.select("chr_strand").unique().to_dicts()
    chrStrand_2_geneInfos = {pair["chr_strand"]: [] for pair in chr_strand_pairs}


    id2geneInfos = {}
    # Trier les gènes par chr_strand puis par position de début
    genes = (
        merged_data.select("gene", "gene_start", "gene_end", "strand", "chr_id", "chr_strand")
        .unique()
        .sort(["chr_strand", "gene_start", "gene_end"])  # Tri par chr_strand d'abord, puis par position
        .to_dicts()
    )
    for row in genes:
        gene_id = row["gene"]
        chr_strand =row["chr_strand"]
        if (gene_id in id2geneInfos):
            print (f"Error: Gene {gene_id} appears multiple times in the GFF file.")
            sys.exit(1)
        geneInfo =GeneInfo(
            chr=row["chr_id"],
            strand=row["strand"],
            gene_id=row["gene"],
            gene_bounds=Bounds(start=row["gene_start"], end=row["gene_end"])
        )
        id2geneInfos[gene_id] = geneInfo
        chrStrand_2_geneInfos[chr_strand].append(geneInfo)
    del genes
       # Option 1: Obtenir tous les mRNAs uniques du DataFrame merged_data
    mRNAs = merged_data.select("mRNA", "gene", "start", "end", "phase", "strand","chr_strand").to_dicts() 
    del merged_data
    mrna_ids = {} 
    # Process the data
    for row in mRNAs:
        # Create MrnaInfo object
        mrna_id=row["mRNA"]
        if mrna_id in mrna_ids:
            print(f"Error: mRNA {mrna_id} appears multiple times in the GFF file.")
            sys.exit(1)
        if row["gene"] not in id2geneInfos:
            print(f"Error: mRNA {mrna_id} has no associated gene.")
            sys.exit(1)
        gene = id2geneInfos[row["gene"]]
        
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
                gene=gene,
                mrna_id=row["mRNA"],
                cds_bounds=cds_coords,
                phase=row["phase"][0] if row["strand"] == "+" else row["phase"][-1],
                mrna_bounds=Bounds(cds_bounds[0].start, cds_bounds[-1].end)
                )
        )

        # Add mRNA ID to the dictionary
        mrna_ids[mrna_id] = 1 
    del mRNAs
    for gene in id2geneInfos.values():
        # define the coding bounds using min and max of the mRNA bounds
        if not gene.mRNAs:
            print(f"Error: Gene {gene.gene_id} has no mRNAs.")
            sys.exit(1)
        gene.coding_bounds = Bounds(
            start=min(mrna.mrna_bounds.start for mrna in gene.mRNAs),
            end=max(mrna.mrna_bounds.end for mrna in gene.mRNAs)
        )
    return chrStrand_2_geneInfos